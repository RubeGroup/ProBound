package modelOptimization;

import proBoundApp.*;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Timestamp;
import java.util.*;

import org.json.*;

import modelComponents.BindingMode;
import modelComponents.BindingModeInteraction;
import modelComponents.CountTable;
import modelComponents.EnrichmentModel;
import modelComponents.ExponentialKineticsModel;
import modelComponents.ModelComponent;
import modelComponents.SELEXModel;

public class LikelihoodOptimizer {
	
	//GOAL: Abstract class containing code:
	// - Running optimization.
    // - Writing trajectory etc.
	// - Knows about what binding modes are updated / optimized.
	// - 
	//
	
	
	// QUESTION: How do we coordinate the fitting strategy?
	// APPROACH:
	//  1) Each model component (binding mode, binding mode interaction, or rho, gamma or h) has:
	//    - boolean included;        //Indicates if this component has been included yet.
	//    - boolean optimized;       //Indicates if the component has been fully optimized yet.
	//    - boolean fit;             //Indicates if the binding mode should be fit.
	// 2) Given this, we use getJSONState() and setStateJSON() to record and update all INCLUDED components.
	// 3) LikelihoodOptimizer creates an ordered indicating when the binding modes should be included.
	//    - Binding modes are ordered in order of ID
	//    - Once the relevant interactions have been included, the interactions are included.
	//    - Potentially we could add one experiment at a time.
	// 4) Each binding component can somehow create a list of possible variations.
	//    - ArrayList<JSONObject> modelComponent.getVariations(JSONObject in) creates a list of JSON objects that encode variations.
	//       - If this returns null when the component is fully optimized. 
	// 5) Add optional fitting strategy:
	//    - Default is to add sequentially/after unlocking. All experiments are added initially.
	//    - In the model specification, add option:
	//           [ {"class":"bindingMode",             "id":0, "initial": true},
	//             {"class":"bindingMode",             "id":1, "initial": false},
	//             {"class":"bindingMode",             "id":2, "initial": false},
	//             {"class":"bindingModeInteractions", "id":0, "initial": false},
	//             {"class":"experiment",              "id":0, "initial": false} ]

	
	//Variables for setting up the numerical optimization object
	////////////////////////////////////////////////////////////
	modelOptimization.Optimizer o;
	String minimizerType = null;
	
	//SGD Settings
	String sgdMethod;
	int sgdBatchsize, sgdMaxIters;
	double sgdLearningRate, sgdConvergence;
	double[] dgdAdaptive;

	//SLBFGS Settings
	int sLGradientBatch, sLHessianBatch, sLMemorySize, sLMaxEpoch, sLHessianPeriod, sLEpochIterations;
	double sLStepSize, sLEpsilon, sLdelta, sLrms_gamma; 
	boolean sLadaptive, sLsvrg;

	// Settings for combined likelihood
	///////////////////////////////////
	double lambdaL2, expBound, likelihoodThreshold;
	int nThreads, nRetries;
	boolean fixedLibrarySize;
	
	//Settings for output
	/////////////////////
	String outputPath, baseName;
	boolean verbose, storeHessian, printPSAM, printTrajectory;
	
	public LikelihoodOptimizer(JSONObject config) {
	
		JSONObject oOptSet = config.getJSONObject("optimizerSetting");
		

		//1. Determines where to write? Makes the CombinedLikelihood  
		if(oOptSet.has("output")) {
			JSONObject oOut = oOptSet.getJSONObject("output");
			verbose         = oOut.getBoolean("verbose");	
			outputPath      = oOut.has("outputPath") ? oOut.getString("outputPath") : null;
			baseName        = oOut.has("baseName")   ? oOut.getString("baseName")   : null;
			storeHessian    = oOut.getBoolean("storeHessian");
			printPSAM       = oOut.getBoolean("printPSAM");
			printTrajectory = oOut.getBoolean("printTrajectory");			
		} else {
			verbose         = false;	
			outputPath      = null;
			baseName        = null;
			storeHessian    = false;
			printPSAM       = false;
			printTrajectory = false;			
		}
		
		//2. Based on config, builds optimizer object. 
		
		minimizerType       = oOptSet.getString("minimizerType");
		lambdaL2            = oOptSet.getDouble("lambdaL2");
		expBound            = oOptSet.getDouble("expBound");
		nThreads            = oOptSet.getInt("nThreads");
		nRetries            = oOptSet.getInt("nRetries");
		likelihoodThreshold = oOptSet.getDouble("likelihoodThreshold");
		
		fixedLibrarySize = oOptSet.getBoolean("fixedLibrarySize"); 

		
		if(		        minimizerType.equals("patternSearch")) {
			if(verbose)
				System.out.println("> Using Pattern Search.");
			JSONObject oSett = oOptSet.getJSONObject("patternSearchSettings");
			o = new modelOptimization.PatternSearch(oSett, verbose);
			
		} else if(		minimizerType.equals("lbfgs")) {
			if(verbose)
				System.out.println("> Using LBFGS.");
			JSONObject oSett = oOptSet.getJSONObject("lbfgsSettings");
			o = new modelOptimization.LBFGS(oSett, verbose);
			
		} else if(		minimizerType.equals("sgd")) {
			if(verbose)
				System.out.println("> Using Stochastic Gradient Descent.");
			JSONObject oSett = oOptSet.getJSONObject("sgdSettings");
			o = new modelOptimization.StochasticGradientDescent(oSett, verbose);
			
		} else {
			System.err.println("ERROR: Invalid optimizer: "+minimizerType);
			System.exit(1);
		}

	}
	
	public void optimizeLikelihood(CombinedLikelihood l) {
		
		
		//1. Setup.
		String fitId = "";
		//1.1 Excludes all components from the model.
		for(int iComp=0; iComp<l.componentList.size(); iComp++ )
			l.componentList.get(iComp).setComponentInclusion(false);
		
		//2. Loop over new components to add.
		for(int iComp=0; iComp<l.fittingOrder.size(); iComp++ ) {
			
			//2.1 Includes the new component.
			ArrayList<ModelComponent> aComp = l.fittingOrder.get(iComp);
			if(aComp.size()==0)
				continue;
			
			ModelComponent oComp =aComp.get(0);
			
			
			String idName = oComp.componentName;
			String li = "";
			for(li=""; li.length()<idName.length(); li += "=");
			System.out.println("");
			System.out.println("===================="+li);
			System.out.println("== Starts fiting "+idName+" ==");
			System.out.println("===================="+li);

			if(oComp instanceof BindingMode || oComp instanceof BindingModeInteraction) {
				//For binding modes and interactions:
				// 1) Include binding mode/interaction, fix all binding modes/interactions, and optimize h
				// 2) Perform initial optimization of the binding mode or interactions (all other components fixed)
				// 3) Explore variations of the initial binding mode/interaction
				// 4) If the added binding mode is frozen, unfreeze it.
				// 5) Perform final fit with all components free. 
				
				if(aComp.size()>1) 
					throw new java.lang.RuntimeException("ERROR: Multi-fitting is not implemented for binding modes.");

				//1.0 Includes binding mode
				oComp.setComponentInclusion(true);
				
				//1.1 Checks if any new experiment can be included.
				oComp.setComponentInclusion(true);
				for(int iExp=0; iExp<l.enrichmentModels.size(); iExp++) {
					EnrichmentModel oEnr = l.enrichmentModels.get(iExp); 
					CountTable oTab      = l.tableModels.get(iExp);
					if(oEnr.includeComponent)
						continue;
					for(int iBM=0; iBM<oEnr.bindingModes.size(); iBM++)
						if(oEnr.bindingModes.get(iBM).includeComponent) {
							oEnr.setComponentInclusion(true);
							oTab.setComponentInclusion(true);
							break;
						}
				}
				
				//2.1.1 Adds component, adjusts the of the component so that it is subdominant.
				oComp.setComponentInclusion(false);
				oComp.activationAdjustment(0.05);
				oComp.setComponentInclusion(true);
								
				//2.1.2 If the new component makes the binding model sequence specific and the enrichment model
				//      has binding saturation, then rescale the activities to make sure a <~ 1
				if(oComp instanceof BindingMode && ((BindingMode) oComp).k > 0) {
					
					//Computes the mean concentration.
					double meanConc = 0.0;
					for(EnrichmentModel oEnr : l.enrichmentModels)
						meanConc += Math.log(oEnr.concentration);
					meanConc /= l.enrichmentModels.size();
					meanConc = Math.exp(meanConc);
					//Loop over binding modes and considers SELEX models with binding saturation.
					for(EnrichmentModel oEnr : l.enrichmentModels) {
						if(oEnr instanceof SELEXModel && ((SELEXModel) oEnr).bindingSaturation) {
							SELEXModel oSELEX = (SELEXModel) oEnr;
							//Checks if any of the old binding modes are sequence specific. 
							boolean enrichmentSpecific = false;
							for(BindingMode oBM : oSELEX.bindingModes) {
								if( oBM == ((BindingMode) oComp) && oBM.includeComponent && oBM.k>0 ) {
									enrichmentSpecific = true;
								}
							}
							
							//If there were no sequence-specific binding mode in the SELEX model, adjust the activities.
							if(enrichmentSpecific) {
								System.out.println(">>> "+oComp.componentName+" is the first sequence-specific componenent in "+
										oSELEX.componentName+" and this enrichment model has bindingSaturation="+
										oSELEX.bindingSaturation+".");
								oSELEX.activationAdjustment(0.05 * oEnr.concentration / meanConc);
							}
						}
					}
				}
				
				//2.3 Fits the count table h-parameter after fixing everything else.
				fitId = "component"+iComp+"-0-h";
				System.out.println("> Optimizing h ("+fitId+").");
				for(int iAll=0; iAll<l.componentList.size(); iAll++) 
					l.componentList.get(iAll).setComponentFiting(false);
				for(int iTab=0; iTab<l.tableModels.size(); iTab++)
					if(l.tableModels.get(iTab).includeComponent)
						l.tableModels.get(iTab).setComponentFiting(true);
				optimizeCurrentModel(l, fitId);

				//2.4 Activates the new component and optimizes.
				if(oComp instanceof BindingMode || oComp instanceof BindingModeInteraction) {
					oComp.setComponentFiting(true);
					fitId = "component"+iComp+"-1-f"+oComp.freezingLevel;
					System.out.println("> Initial optimization ("+fitId+").");
					optimizeCurrentModel(l, fitId);
				}

				//2.5 Explores variations of the new model component.
				double logL_best                       = l.value;
				JSONObject bestState                   = l.getJSONState();
				String bestKey                         = null;

				Hashtable<String,JSONObject> variation = oComp.getVariations(bestState, l.componentList);
				int nVariation                         = 0;


				while(variation!=null) {

					System.out.println("Suggested variations:");
					for(String k : variation.keySet())
						System.out.println("key="+k+", description = "+oComp.variationDescription.get(k));


					String impovedName          = null; 

					//2.5.1 Loops over the suggested variations.
					for(String key : variation.keySet()) {
						//2.6.1.1 Sets the state to the current variation.
						l.setStateJSON(variation.get(key));

						//2.5.1.2 Optimizes
						fitId = "component"+iComp+"-2-variation"+nVariation;
						System.out.println("> Optimizing variation \""+oComp.variationDescription.get(key)+"\" ("+fitId+").");
						optimizeCurrentModel(l, fitId);
						nVariation++;

						//2.5.1.3 Saves the model if the likelihood improved. Otherwise removes the trajectory file.
						if(l.value < logL_best - likelihoodThreshold) {
							System.out.println("  The Likelihood DID improve.");
							//2.5.1.3.1 Removes the old (improved) trajectory file
							if(impovedName!=null)
								removeTrajectory(impovedName);

							//2.5.1.3.2 Saves the current state.
							logL_best   = l.value;
							bestState   = l.getJSONState();
							bestKey     = key;
							impovedName = fitId;
						} else {
							System.out.println("  The Likelihood DID NOT improve. Discarding fit "+fitId+".");
							//2.5.1.3.3 Removes the trajectory file if the likelihood didn't improve
							removeTrajectory(fitId);
						}
					}

					//2.5.2 Set the likelihood to the best state.
					l.setStateJSON(bestState);
					oComp.indicateModelChoice(bestKey);

					//2.5.3 Gets new variations...
					variation = oComp.getVariations(bestState, l.componentList);

				}

				//2.6 Unfreezes the new model component
				while(oComp.freezingLevel>0) {
					fitId = "component"+iComp+"-3-f"+(oComp.freezingLevel-1);
					System.out.println("> Unfreezing component ("+fitId+").");
					oComp.setFreezingLevel(oComp.freezingLevel-1);
					optimizeCurrentModel(l, fitId);
				}

			} else if (oComp instanceof SELEXModel){
				
				// Strategy for optimizing the SELEX model:
				// 1) Save the current model
				// 2) Get variations of the SELEX model component.
				// 3) Optimize h for the count table attached to the SELEX model.
				// 4) Let all components free, fit
				// 5) Compares the likelihood with the initial model and only keep if it improved.
				
				// 2.1 Saves the last state.
				double logL_best                       = l.value;
				JSONObject bestState                   = l.getJSONState();
				String bestKey                         = null;
				String impovedName                     = null; 
				String impovedNameH                    = null; 

				// 2.2 Gets (one) variation of the SELEX model
				Hashtable<String,JSONObject> variation = oComp.getVariations(bestState, l.componentList);
				System.out.println(variation);
				int nVariations                        = variation.size();
				int iVariation                         = 0;
				if(nVariations==0)
					continue;

				if(nVariations>1)
					System.err.println("WARNING: Expected at most one variation of '"+oComp.componentName+"' but "+variation.size()+" were generated. Proceeding by testing all variations and keeping the best. ");
				
				for(String key : variation.keySet()) { 
					l.setStateJSON(variation.get(key));
					//The activity is shifted so that the expected activity should is 5% binding.
					oComp.activationAdjustment(0.05);
					
					// 2.3 Fits the count table h after fixing all other components. 
					String fitIdH = "component"+iComp+"-1-variation-"+iVariation+"-0-h";
					System.out.println("> Optimizing h ("+fitIdH+").");
					for(int iAll=0; iAll<l.componentList.size(); iAll++) 
						l.componentList.get(iAll).setComponentFiting(false);
					for(int iTab=0; iTab<l.tableModels.size(); iTab++)
						if(l.tableModels.get(iTab).includeComponent)
							l.tableModels.get(iTab).setComponentFiting(true);
					optimizeCurrentModel(l, fitIdH);
					
					// 2.4 Fits all the full model 
					fitId = "component"+iComp+"-1-variation-"+iVariation+"-1-all";
					System.out.println("> Optimizing the full model ("+fitId+").");
					for(int iAll=0; iAll<l.componentList.size(); iAll++)
						if(l.componentList.get(iAll).includeComponent)
							l.componentList.get(iAll).setComponentFiting(true);
					optimizeCurrentModel(l, fitId);
					
					// 2.5) Compares the likelihood with the initial model and only keep if it improved.
					if(l.value < logL_best) {
						System.out.println("  The Likelihood DID improve.");
						//2.5.1.3.1 Removes the old (improved) trajectory file
						if(impovedName!=null){
							removeTrajectory(impovedName);
							removeTrajectory(impovedNameH);
						}

						//2.5.1.3.2 Saves the current state.
						logL_best    = l.value;
						bestState    = l.getJSONState();
						bestKey      = key;
						impovedName  = fitId;
						impovedNameH = fitIdH;

					} else {
						System.out.println("  The Likelihood DID NOT improve. Discarding fit "+fitId+".");
						//2.5.1.3.3 Removes the trajectory file if the likelihood didn't improve
						removeTrajectory(fitId);
						removeTrajectory(fitIdH);
					}
				}

				//2.5.2 Set the likelihood to the best state.
				l.setStateJSON(bestState);
				oComp.indicateModelChoice(bestKey);
				
			} else if (oComp instanceof ExponentialKineticsModel) {
				//Checks so all enrichment models are of the same type
				for(ModelComponent compI : aComp)
					if(!(compI instanceof ExponentialKineticsModel)) 
						System.err.println("WARNING: All component in this fitting batch must be ExponentialKineticsModel, but componentName= "+oComp.componentName);
					
				//Identifies what freezing level to start with.
				int maxLevel = 0;
				for(ModelComponent compI : aComp)
					maxLevel = Math.max(maxLevel, compI.freezingLevel);

				while(maxLevel>0) {
					fitId = "component"+iComp+"-3-f"+(oComp.freezingLevel-1);
					System.out.println("> Unfreezing components ("+fitId+").");

					//Unfreezes the most frozen components
					for(ModelComponent compI : aComp) {
						if(compI.freezingLevel == maxLevel) {
							System.out.println("   Component: "+oComp.componentName+".");

							compI.setFreezingLevel(compI.freezingLevel-1);
						}
					}

					//Optimizes current batch
					optimizeCurrentModel(l, fitId);

					//Updates the maximum freezing value
					maxLevel=0;
					for(ModelComponent compI : aComp) { 
						maxLevel = Math.max(maxLevel, compI.freezingLevel);
					}
					

				}
				
				
			} else {
				System.err.println("WARNING: No optimization of the model component '"+oComp.componentName+"' is implemented.");
			}
			
			//2.7 Optimizes all included components
			fitId = "component"+iComp+"-4-all";
			System.out.println("> Optimizing the full model ("+fitId+").");
			for(int iAll=0; iAll<l.componentList.size(); iAll++)
				if(l.componentList.get(iAll).includeComponent)
					l.componentList.get(iAll).setComponentFiting(true);
			optimizeCurrentModel(l, fitId);
			
		}
		
	}
	
	private void removeTrajectory(String name) {
		if(printTrajectory && outputPath!=null && baseName!=null) {
			String trajectoryFile = outputPath+"/"+baseName+".trajectory."+name+".csv";
			try {
				Files.delete(Paths.get(trajectoryFile));
			} catch (Exception e) {
				System.out.println("No trajectory file to delete..");
			}
		}
	}
	
	public void optimizeCurrentModel(CombinedLikelihood l, String optName) /*throws Exception*/ {
		
		System.out.println(">>  Starting new optimization: "+optName+". ("+(new Timestamp(System.currentTimeMillis())).toString()+").");
		l.updatePacking();
		
		System.out.println(">>> Packing before optimization");

		///////////////////////////
		// Computes gradient before 
		l.updatePacking();
		//l.lossFunction_updateGradient();
		//l.bindingModes.get(1).setFreezingLevel(0);
		//System.out.println("Packing:     "+l.packing.getJSONObject("packing").toString());
		//System.out.println("");
/*		System.out.println("================= GRADIENT BEFORE OPTIMIZATION =================");
		l.lossFunction_updateGradient();
		double[] vGrad  = l.gradient;
		double[] vFGrad = l.gradientFiniteDifferences(0.001);
		System.out.println("Gradient:         "+Misc.formatVectorE_d(l.gradient, ",", "{", "}", 8));
		System.out.println("FD Gradient:      "+Misc.formatVectorE_d(vFGrad, ",", "{", "}", 8));
		System.out.println("|Grad - FD Grad|: "+Array.norm(Array.subtract(vGrad, vFGrad)));
		System.out.println("Gradient norm:    "+Array.norm(vGrad));
		System.out.println("================================================================");
		System.out.println("");*/
		///////////////////////////
		
		System.out.println("Value and gradient before optimization:");
		System.out.println("=======================================");
		l.lossFunction_updateGradient();
		System.out.println("value         = "+l.value);
		System.out.println("gradient      = "+Misc.formatVector_d(l.gradient));
		System.out.println("gradient norm = "+Array.norm(l.gradient));
		
		
		
		if(printTrajectory && outputPath!=null && baseName!=null) 
			l.startNewTrajectoryFile(outputPath+"/"+baseName+".trajectory."+optName+".csv", l.packing);
		
		for(int iTry=0; iTry<nRetries; iTry++) {
			try {
				l.lossFunction_reportPosition();
				o.optimize(l);
//				break;
			} catch (Optimizer.Optimizer_MaxIteractionException e) {
				System.out.println(">>> Maximum iteration count reached.");
				
				/*System.out.println("================= CHECKING GRADIENT USING FD ==================");
				l.lossFunction_updateGradient();
				double[] vGrad  = l.gradient;
				double[] vFGrad = l.gradientFiniteDifferences(0.001);
				System.out.println("Gradient:         "+Misc.formatVectorE_d(l.gradient, ",", "{", "}", 8));
				System.out.println("FD Gradient:      "+Misc.formatVectorE_d(vFGrad, ",", "{", "}", 8));
				System.out.println("|Grad - FD Grad|: "+Array.norm(Array.subtract(vGrad, vFGrad)));
				System.out.println("Gradient norm:    "+Array.norm(vGrad));
				System.out.println("================================================================");*/
				
				
//				break;
				
			} catch (Exception e) {
	        	e.printStackTrace();
	        	if(l.lastReported != null)
	        		l.lossFunction_setParameters(l.lastReported);
	        	System.out.println(">>>> Exception caught. Parameters reverted.");
	        	System.out.println("> Parameter: "+Misc.formatVectorE_d(l.parameters));
	        	System.out.println("> Gradient:  "+Misc.formatVectorE_d(l.gradient));
	        	if(iTry<nRetries) {
	        		System.out.println(">>>> Re-trying ("+(iTry+1)+"/"+nRetries+").");
	        		continue;
	        	} else {
	    			System.out.println(">>>> Optimization failed. Breaks.");
	        		break;
	        	}
			}
			
			break;
		}
		
		System.out.println("After: gradient norm = "+Array.norm(l.gradient));
		
		printModelToJSONList(l, optName);
		
		System.out.println(">>> Parameters after optimization");
		JSONModel.printJSONState(l.getJSONState(), "coefficients");
		System.out.println("");
		//System.out.println("gradient = "+Misc.formatVector(l.gradient));
				
		if(printTrajectory)
			l.stopTrajectoryOutput();
	}
	
	public void printModelToJSONList(CombinedLikelihood l, String optName) {

		JSONObject currentState = l.getJSONState();	
		//Adds the optimizerSetting:
		JSONObject optimizerSetting = new JSONObject();
		
		//Saves the optimizer settings.
		currentState.put("optimizerSetting",        optimizerSetting);
		optimizerSetting.put("lambdaL2",            lambdaL2);
		optimizerSetting.put("expBound",            expBound);
		optimizerSetting.put("nThreads",            nThreads);
		optimizerSetting.put("nRetries",            nRetries);
		optimizerSetting.put("likelihoodThreshold", likelihoodThreshold);
		optimizerSetting.put("fixedLibrarySize",    fixedLibrarySize);
		
		if(o!=null)
			optimizerSetting.put(minimizerType+"Settings", o.saveToJSON_settings());
		
		JSONObject output = new JSONObject();
		optimizerSetting.put("output", output);
		if(outputPath!=null)
			output.put("outputPath", outputPath);
		if(baseName!=null)
			output.put("baseName", baseName);
		
		//Adds the metadata.
		currentState.put("metadata", l.getMetadataJSON(null, optName));

		//Saves the modelComponent metadata.
		for(int iC=0; iC<l.componentList.size(); iC++)
			l.componentList.get(iC).saveToJSON_settings_metadata(currentState);
		
		if(outputPath!=null && baseName!=null)
			writeCompactJSONModel(currentState,outputPath+"/"+baseName+".models.json", true);
	}
	
	//Writes the current model as a single-line JSON object. 
	public static void writeCompactJSONModel(JSONObject model, String outPath, boolean append) {
		PrintStream original	= System.out;
		
		//Changes output stream.
		if(outPath != null){
			try {
				System.setOut(new PrintStream(new FileOutputStream(outPath, append)));
				
			} catch (FileNotFoundException e) {
				System.out.println("Cannot create trajectory file at this "
						+ "location: "+outPath);
				e.printStackTrace();
			}
		
		}
		
		System.out.println(model.toString());
		
		//Returns to original stream.
		System.setOut(original);
	}
	
}

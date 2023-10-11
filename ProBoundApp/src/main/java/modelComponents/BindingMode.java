package modelComponents;

import java.util.*;
import org.json.*;

import base.Array;
import sequenceTools.LongSequence;
import proBoundApp.Misc;
import sequenceTools.SlidingWindow;

public class BindingMode extends ModelComponent  {
	
	//Model settings.
	public int k, flankLength, dInt, dIntMax;
	public String symmetryString;
	//BindingModeSymmetry symmetry;
	public boolean usePositionBias, singleStrand;
	
	//Variables for fixing shifts and
	private int currentShift;
	
	//Indicators for determining what sub-components to fit.
	public boolean fitMono, fitDi, fitActivity, fitPositionBias;
	
	//Fitting strategy and constraints
	public int wPositionBiasBin, minSize, maxSize, maxFlankLength;
	public boolean fitK, fitFlankLength, shiftMotif, heuristicShift, optimizeSizeHeuristic;
	public boolean roundSpecificActivity, experimentSpecificPositionBias, experimentSpecificActivity;
	public boolean fitLogActivity;
	JSONArray fittingStages;
	public int iFittingStage;
	public double informationThreshold;
	
	//Binding specificity parameters
	public double[] monoBetas;
	public ArrayList<double[]> diBetas;
	
	//Variables for holding 
	public ArrayList<String> modifications;                 //List of modifications
	public ArrayList<Character[]> monoMod;           //List of modified mononucleotides
	public ArrayList<ArrayList<Character[][]>> diMod;//List of modified dinucleotides
	public ArrayList<int[]> monoBetaShiftIndices; //Indices and values of the components in monoBetas and diBetas that should be shifted.
	public ArrayList<ArrayList<int[]>> diBetaShiftIndices;
	public ArrayList<double[]> monoBetaShiftValues;   
	public ArrayList<ArrayList<double[]>> diBetaShiftValues;
	
	//Experiment-specific parameters
	public ArrayList<CountTable> countTables;
	public ArrayList<Integer> maxFrames;
	public ArrayList<ArrayList<double[]>> positionBiasBetas, positionBiasAlphas;
	public ArrayList<double[]> activityBetas, activityAlphas;
	
	public boolean swIncludeDi; //Indicates if the optimization removing dinucleotides is used.
	
	LongSequence.SequenceClass sc;
	public String letterComplement, letterOrder;
	public int nMono, nDi;
	
	public BindingMode(JSONObject config, int iBMIn, LongSequence.SequenceClass scIn, String letterComplementIn, String letterOrderIn) {
		super("bindingModes");

		
		iComp           = iBMIn;		
		componentName   = "Binding mode "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");
		
		//Initially we have no experiments.
		letterComplement     = letterComplementIn;
		letterOrder          = letterOrderIn;
		nMono                = letterOrder.length(); 
		nDi                  = nMono*nMono;
		
		countTables           = new ArrayList<CountTable>();
		sc                    = scIn;
		readFromJSON_settings(config);
		iFittingStage         = 0;
		readFromJSON_constraints(config);
		seed_component(config);
		
		maxFreezeLevel = (k>0&&usePositionBias ? 1 : 0) + (dInt>0 ? 1 : 0);
		setFreezingLevel(maxFreezeLevel);
		
		//Initiates system for keeping track of variations.
		currentShift         = 0;
		variationName        = k+";"+flankLength+";"+currentShift;
		variationsOptimized  = false;
		
		diBetaShiftValues    = new ArrayList<ArrayList<double[]>>();
		monoBetaShiftValues  = new ArrayList<double[]>();   
		
		
	}
	
	private void setFittingStage(int iStage) {
		iFittingStage         = iStage;
		JSONObject oCurr      = fittingStages.getJSONObject(iFittingStage);
		fitK                  = oCurr.has("optimizeSize")                ? oCurr.getBoolean("optimizeSize")                        : false;
		fitFlankLength        = oCurr.has("optimizeFlankLength")         ? oCurr.getBoolean("optimizeFlankLength")                 : false;
		optimizeSizeHeuristic = oCurr.has("optimizeSizeHeuristic")       ? oCurr.getBoolean("optimizeSizeHeuristic")               : false;
		shiftMotif            = oCurr.has("optimizeMotifShift")          ? oCurr.getBoolean("optimizeMotifShift")                  : false;
		heuristicShift        = oCurr.has("optimizeMotifShiftHeuristic") ? oCurr.getBoolean("optimizeMotifShiftHeuristic")         : false;

	}
	
	public void allocateParameters() {
		
		monoBetas           = zero_d(nMono*k);
		diBetas             = new ArrayList<double[]>();
		for(int iDi=0; iDi<dInt; iDi++)
			diBetas.add(      zero_d(nDi*(k-(iDi+1))));
		
		maxFrames           = new ArrayList<Integer>();
		activityBetas       = new ArrayList<double[]>();
		activityAlphas      = fitLogActivity ? null : new ArrayList<double[]>();
		if(usePositionBias && k>0)
			positionBiasBetas   = new ArrayList<ArrayList<double[]>>();
		else
			positionBiasBetas   = null;
		
		//Experiment-specific parameters
		for(int iExp=0; iExp<countTables.size(); iExp++) {
			
			maxFrames.add(             2*(countTables.get(iExp).l + 2*flankLength - k + 1) );
			
			//Creates activity variables
			if(fitLogActivity) {
				activityBetas.add(zero_d(countTables.get(iExp).nColumns));
			} else {
				activityBetas.add(null);
				activityAlphas.add(zero_d(countTables.get(iExp).nColumns));
			}
			
			//Creates position bias variables. 
			if(usePositionBias && k>0)
				positionBiasBetas.add( zero_Ad(2, maxFrames.get(iExp)/2) );

		}
		
		//Allocates values for the values for chemical modifications.
		monoBetaShiftValues = new ArrayList<double[]>();
		for(int i=0; i<monoBetaShiftIndices.size(); i++)
			monoBetaShiftValues.add(new double[monoBetaShiftIndices.get(i).length]);
		diBetaShiftValues   = new ArrayList<ArrayList<double[]>>();
		for(int i=0; i<diBetaShiftIndices.size(); i++) {
			ArrayList<int[]>    tempIn  = diBetaShiftIndices.get(i);
			ArrayList<double[]> tempOut = new ArrayList<double[]>(); 
			for(int j=0; j<tempIn.size(); j++)
				tempOut.add(new double[tempIn.get(j).length]);
			diBetaShiftValues.add(tempOut);
		}
		
		updateAlphas();
	}
	

	@Override
	public void setComponentFiting(boolean fit) {
		
		if(fit) {
			if(includeComponent)
				fitComponent = true;
			else
				fitComponent = false;
		} else {
			fitComponent = false;
		}
		
		if(fitComponent) {
			
			//freezingLevels:
			// 1) The position bias is frozen
			// 2) The position bias and dinucleotide interactions are both frozen.
			
			fitMono             = true;
			
			int fDiUnfreeze = k>0&&usePositionBias ? 2 : 1;
			fitDi = dInt>0 && freezingLevel<fDiUnfreeze ? true : false;
			
			if(freezingLevel<1)
				fitPositionBias = usePositionBias && k>0 && freezingLevel<1 ? true : false;

			fitActivity         = true;

		} else {
			fitMono         = false;
			fitDi           = false;
			fitActivity     = false;
			fitPositionBias = false;
		}
	}
	
	public void updateAlphas() {
		//Computes alphas if the fitted parameters is beta=log(alpha)
		if( fitLogActivity )
			activityAlphas     = exp_Ad(activityBetas);
		
		//Updates position bias alphas 
		if(usePositionBias && k>0) 
			positionBiasAlphas = exp_AAd(positionBiasBetas);
		else
			positionBiasAlphas = null;
		
	}
	
	public SlidingWindow getSlidingWindow(int iExp, Set<String> enrMod) {
		
		if(monoBetas.length > 0) {

			//Packs the scoring matrix into a format that the SlidingWindow class can read. 
			ArrayList<ArrayList<double[]>> aal = new ArrayList<ArrayList<double[]>>();
			aal.add(new ArrayList<double[]>());

			//Global betas.
			double[] tempMono = Array.clone(monoBetas);
			//Adds modification shifts to compute the effective betas
			for(int iMod=0; iMod<modifications.size(); iMod++) {
				if(enrMod.contains(modifications.get(iMod))) {
					for(int iParam=0; iParam<monoBetaShiftValues.get(iMod).length; iParam++) {
						tempMono[  monoBetaShiftIndices.get(iMod)[iParam]] 
								+= monoBetaShiftValues.get( iMod)[iParam];
					}
				}
			}
			
			//Saves the effective betas
			aal.get(0).add(tempMono);

			//The dinuceotide coefficients should only be included if 1) fitDi=true
			//or 2) at any diBeta is not zero.
			swIncludeDi = fitDi;
			for(int id=0; id<diBetas.size(); id++) {
				double[] diBetasI   = diBetas.get(id);
				for(int ix=0; ix<diBetasI.length; ix++)
					if(diBetasI[ix]!=0) {
						swIncludeDi = true;
						break;
					}
				if(swIncludeDi)
					break;
			}

			if( dInt>0 && swIncludeDi) {
				aal.add(new ArrayList<double[]>());
				for(int id=0; id< dInt; id++) {
					//Global Betas
					double[] tempDi = Array.clone(diBetas.get(id));
					//Adds modification shifts to compute the effective betas
					for(int iMod=0; iMod<modifications.size(); iMod++) {
						if(enrMod.contains(modifications.get(iMod)))
							for(int iParam=0; iParam<diBetaShiftValues.get(iMod).get(id).length; iParam++) {
								tempDi[  diBetaShiftIndices.get(iMod).get(id)[iParam]] 
									  += diBetaShiftValues.get( iMod).get(id)[iParam];
							}
					}
					//Saves the effetive betas
					aal.get(1).add(tempDi);
				}
			}

			//Extracts the flanks
			String lf = countTables.get(iExp).fullTable.leftFlank;
			String rf = countTables.get(iExp).fullTable.rightFlank;

			//Creates and adds new sliding window.
			return new SlidingWindow(
					lf.substring(lf.length()-flankLength,lf.length()), 
					rf.substring(0,flankLength)
					, sc, aal, letterOrder);
		} else {
			return null;
		}
	}

	public void addCountTable(CountTable exp) {
		countTables.add(exp);
		allocateParameters();
	}
	
	public double computeExpectedAlpha() {
		
		//Mean alpha
		double meanAlpha = 1;
		
		//Product over expected mononucleotide alphas
		for(int q=0; q<k;  q++)
			meanAlpha *=     Array.sum(Array.exp(Arrays.copyOfRange(monoBetas, nMono*q, nMono*(q+1) )))/ nMono;

		//Product over expected dinucleotide alphas
		for(int iD=0; iD<diBetas.size(); iD++)
			for(int q=0; q<k-(iD+1);q++)
				meanAlpha *= Array.sum(Array.exp(Arrays.copyOfRange(diBetas.get(iD), nDi*q, nDi*(q+1))))/nDi;
		
		return meanAlpha;
	}
	
	//Computes the value of the Dirichlet regularization
	public double getDirichletValue(double pseudocount) {

		double out = 0.0;
		int nMono  = letterOrder.length();

		//Computes the weight
		double weight = 0, nIncExp = 0;
		for(CountTable ct: countTables) {
			if(ct.includeComponent) {
				weight += pseudocount / ct.nReads;
				nIncExp += 1;
			}
		}
		weight /= nIncExp;
		
		//Mononucleotide contribution
		if(fitMono) {
			for(int x=0; x<k; x++) {
				double alphaSum = 0;
				for(int iMono=0; iMono<nMono; iMono++) {
					double beta = monoBetas[x*nMono + iMono]; 
					out        -= weight*beta;
					alphaSum   += Math.exp(beta);
				}
				out += weight * nMono * Math.log(alphaSum);
			}
		}
		
		//Dinucleotide contribution
		if(fitDi) {
			int nDi = nMono * nMono;
			for(int iD=0; iD<dIntMax; iD++) {
				for(int x=0; x<k-iD-1; x++) {
					double alphaSum = 0;
					for(int iDi=0; iDi<nDi;iDi++) {
						double beta = diBetas.get(iD)[x*nDi + iDi];
						out        -= weight*beta;
						alphaSum   += Math.exp(beta);
					}
					out += weight * nDi * Math.log(alphaSum);
				}
			}
		}
		
		return out;
	}
	
	public void addDirichletGradient(JSONObject oGradBM, double pseudocount) {
		
		int nMono  = letterOrder.length();
		
		//Computes the weight
		double weight = 0, nIncExp = 0;
		for(CountTable ct: countTables) {
			if(ct.includeComponent) {
				weight += pseudocount / ct.nReads;
				nIncExp += 1;
			}
		}
		weight /= nIncExp;
		
		//Mononucleotide contribution
		if(fitMono) {
			JSONArray aMono  = oGradBM.getJSONArray("mononucleotide");
			double[] expMono = new double[nMono];
			
			for(int x=0; x<k; x++) {
				double expSum = 0;
				for(int iMono=0; iMono<nMono; iMono++) {
					double expBeta = Math.exp(monoBetas[x*nMono + iMono]); 
					expMono[iMono] = expBeta;
					expSum        += expBeta;
				}
				
				for(int iMono=0; iMono<nMono; iMono++) {
					int i = x*nMono + iMono;
					aMono.put(i, aMono.getDouble(i) + weight * (nMono * (expMono[iMono] / expSum) - 1));
				}
			}
		}
		
		//Dinucleotide contribution
		if(fitDi) {
			int nDi        = nMono * nMono;
			double[] expDi = new double[nDi];
			
			for(int iD=0; iD<dInt; iD++) {
				JSONArray aDi  = oGradBM.getJSONArray("dinucleotide").getJSONArray(iD);
				for(int x=0; x<k-iD-1; x++) {

					double expSum = 0;
					for(int iDi=0; iDi<nDi; iDi++) {
						double expBeta = Math.exp(diBetas.get(iD)[x*nDi + iDi]); 
						expDi[iDi]     = expBeta;
						expSum        += expBeta;
					}
					
					for(int iDi=0; iDi<nDi; iDi++) {
						int i = x*nDi + iDi;
						aDi.put(i, aDi.getDouble(i) + weight * (nDi * (expDi[iDi] / expSum) - 1));
					}
				}
			}
		}
	}
	
	@Override
	public void activationAdjustment(double weight) {
		
		//Computes the mean sequence score of this binding mode
		double alphaSeq = computeExpectedAlpha();
	

		for(int iExp=0; iExp<countTables.size(); iExp++) {
			CountTable ct = countTables.get(iExp);
			
			//Computes the expected partition function  
			double[] expectedPartitionFunction = ct.computeExpectedPartitionFunction();
			
			//If  partition function=0 (which it is for the first binding mode), set the partition function to 1. 
			for(int i=0; i<expectedPartitionFunction.length; i++)
				if(expectedPartitionFunction[i]==0)
					expectedPartitionFunction[i] = 1;
			
			
			//Computes the (possibly position-bias weighted) number of binding of frames. 
			double alphaPBSum = (usePositionBias && k>0) ? tr_Ad(positionBiasAlphas.get(iExp)) : maxFrames.get(iExp);

			//Computes the appropriate activity rescaling factor.
			double[] rescalingFactor =  Array.scalarMultiply(expectedPartitionFunction, weight / (alphaSeq * alphaPBSum * ct.enr.concentration));
			
			//Shifts the activity
			//NOTE: If we have roundSpecificActivity=false or experimentSpecificActivity-false, we need to make sure
			//that the adjustment keeps the consistent with these settings. Note that the naive rescaling factor should 
			//be constant if these settings are used.   
			for(int iCol=0; iCol<activityBetas.get(iExp).length; iCol++) {
				if(rescalingFactor[iCol]==0) //OLD: The partition function will be zero before the first binding mode is added (should not happen).  
					continue;
				if(fitLogActivity)
					activityBetas.get(iExp)[iCol] += Math.log(rescalingFactor[iCol]);
				else
					activityAlphas.get(iExp)[iCol] *= rescalingFactor[iCol];
						
			}
		}


		if(fitLogActivity)
			updateAlphas();

	}
	
	@Override
	public Hashtable<String,JSONObject> getVariations(JSONObject currentModel, ArrayList<ModelComponent> componentList) {
		
		if(variationsOptimized)
			return null;
		
		Hashtable<String,JSONObject> newVariations = new Hashtable<String,JSONObject>();
		String testLabel;
		
		//1. Determines if the binding mode allows for expanded footprint or increased flank length
		
		//Determines if it is possible to expand the footprint of the binding mode. 
		boolean expandFootprintPossible =  ( maxSize<0 || (maxSize>0 && (k+1)<=maxSize) );
		for(int iExp=0; iExp<countTables.size(); iExp++) {
			CountTable ct = countTables.get(iExp);
			int maxSeqLength;
			if(maxFlankLength < 0) 
				maxSeqLength = ct.l + 2 * Math.min(ct.leftFlank.length(), ct.rightFlank.length());
			else 
				maxSeqLength = ct.l + 2 * maxFlankLength;
			if((k+1) > maxSeqLength)
				expandFootprintPossible = false;
		}
		
		//Determines if expansion of flanks should be tried
		boolean increasedFlankLengthPossible =  
				( maxFlankLength<0 || (maxFlankLength>0 && flankLength<maxFlankLength) )
				&& flankLength < k-1;
		for(int iExp=0; iExp<countTables.size(); iExp++) 
			if(flankLength >= Math.min(
					countTables.get(iExp).rightFlank.length(),
					countTables.get(iExp).leftFlank.length()))
				increasedFlankLengthPossible = false;
		
		//2. Generates variations  
		while(true) {
			
			//2.1 Generates new variations in the current fitting stage.
			//Initial model shift
			/////////////////
			testLabel = k+";"+flankLength+";"+(currentShift);
			if(!variationDescription.containsKey(testLabel)) {
				variationDescription.put(testLabel, "Initial model.");
				newVariations.put(testLabel, currentModel);
			}
			if(newVariations.size()>0)
				return newVariations;

			//Heuristic shift
			/////////////////
			if(heuristicShift && k>1) {
				double centerPos = computeInformationWeightedMean();
				int shift        = - (int) Math.round(centerPos - 0.5*(k-1));
				testLabel        = k+";"+flankLength+";"+(currentShift+shift);
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Heuristic shift ("+shift+").");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation,  shift, -shift, 0, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				if(newVariations.size()>0)
					return newVariations;
			}

			//Increases flank length
			////////////////////////
			if(fitFlankLength && increasedFlankLengthPossible) {
				testLabel = k+";"+(flankLength+1)+";"+currentShift;
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Increases flank length.");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation,  0,  0, 1, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				if(newVariations.size()>0)
					return newVariations;
			}

			//Shifting to the left and right
			////////////////////////////////
			if(shiftMotif && k>1) {
				testLabel=k+";"+flankLength+";"+(currentShift+1);
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Shifts to the right.");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation,  1, -1, 0, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				testLabel=k+";"+flankLength+";"+(currentShift-1);
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Shifts to the left.");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation, -1,  1, 0, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				if(newVariations.size()>0)
					return newVariations;
			} 
			
			//Tests extending the binding mode to the left and right.
			/////////////////////////////////////////////////////////
			if(fitK && expandFootprintPossible) {
				testLabel=(k+1)+";"+flankLength+";"+(currentShift+1);
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Expands to the left.");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation,  1,  0, 0, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				testLabel=(k+1)+";"+flankLength+";"+(currentShift);
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Expands to the right.");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation,  0,  1, 0, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				if(newVariations.size()>0)
					return newVariations;				
			} 
			
			//Tests both extending and increasing the flank length
			/////////////////////////////////////////////////////
			if(fitK && fitFlankLength && expandFootprintPossible && increasedFlankLengthPossible) {
				testLabel=(k+1)+";"+(flankLength+1)+";"+(currentShift+1);
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Expands to the left and increases the flank length.");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation,  1,  0, 1, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				testLabel=(k+1)+";"+(flankLength+1)+";"+(currentShift);
				if(!variationDescription.containsKey(testLabel)) {
					variationDescription.put(testLabel, "Expands to the right and increases the flank length.");
					JSONObject modelVariation = clone_JSON_O(currentModel);
					modifyBindingMode(modelVariation,  0,  1, 1, dIntMax, componentList);
					newVariations.put(testLabel, modelVariation);
				}
				if(newVariations.size()>0)
					return newVariations;
			} 
			
			//Tests both extending and increasing the flank length
			/////////////////////////////////////////////////////
			if(optimizeSizeHeuristic && expandFootprintPossible) {
				//Determines if the optimizer should try to expand to the left and the right
				int deltaLeft = 0, deltaRight=0;
				
				//The binding mode is only expanded to the left/right if the left/right-most 2 bases
				//have information above a certain threshold.
				System.out.println("informationThreshold = "+informationThreshold);
				if(computeBaseInformation(0)>informationThreshold)
					deltaLeft  = 1;
				else if(k>1 && computeBaseInformation(1)>informationThreshold)
					deltaLeft  = 1;
				if(computeBaseInformation(k-1)>informationThreshold)
					deltaRight = 1;
				else if(k>1 && computeBaseInformation(k-2)>informationThreshold)
					deltaRight = 1;
				//By default, flankLength is incremented by one.
				int deltaFlank = increasedFlankLengthPossible ? 1 : 0;
				
				//Determines if it is possible to expand the binding frame by two bases.
				boolean doubleExpansionPossible =  ( maxSize<0 || (maxSize>0 && (k+2)<=maxSize) );
				
				System.out.println(">>> Information used for optimizeSizeHeuristic: "+
						"I[0] = "+computeBaseInformation(0) + 
						", I[1] = "+computeBaseInformation(1)+
						", I[k-2] = "+computeBaseInformation(k-2)+
						", I[k-1] = "+computeBaseInformation(k-1));
				
				for(int iExp=0; iExp<countTables.size(); iExp++) {
					CountTable ct = countTables.get(iExp);
					int maxSeqLength;
					if(maxFlankLength < 0) 
						maxSeqLength = ct.l + 2 * Math.min(ct.leftFlank.length(), ct.rightFlank.length());
					else 
						maxSeqLength = ct.l + 2 * maxFlankLength;
					if((k+2) > maxSeqLength)
						doubleExpansionPossible = false;
				}
				
				//Generates variations of the model.
				if((deltaLeft+deltaRight)==2 && !doubleExpansionPossible ) {
					//In this case 1) the information suggests that we should expand to both to the left and to the 
					//right, and 2) the model settings allow for at most one expansion. We therefore generate
					//two variations, one expanding to the left and one to the right.
					testLabel=(k+1)+";"+(flankLength+1)+";"+(currentShift+1);
					if(!variationDescription.containsKey(testLabel)) {
						variationDescription.put(testLabel, "Expands to the left and increases the flank length.");
						JSONObject modelVariation = clone_JSON_O(currentModel);
						modifyBindingMode(modelVariation,  1,  0, 1, dIntMax, componentList);
						newVariations.put(testLabel, modelVariation);
					}
					testLabel=(k+1)+";"+(flankLength+1)+";"+(currentShift);
					if(!variationDescription.containsKey(testLabel)) {
						variationDescription.put(testLabel, "Expands to the right and increases the flank length.");
						JSONObject modelVariation = clone_JSON_O(currentModel);
						modifyBindingMode(modelVariation,  0,  1, 1, dIntMax, componentList);
						newVariations.put(testLabel, modelVariation);
					}
					
				} else if((deltaLeft+deltaRight)==1 || (doubleExpansionPossible && (deltaLeft+deltaRight)==2)) {
					//In this case, we expand to the left and/or to the right based on the information analysis above.
					testLabel=(k+deltaLeft+deltaRight)+";"+(flankLength+deltaFlank)+";"+(currentShift+deltaLeft);
					if(!variationDescription.containsKey(testLabel)) {
						variationDescription.put(testLabel, "Expands: deltaLeft = "+deltaLeft+", deltaRight = "+deltaRight+", deltaFlank = "+deltaFlank+".");
						JSONObject modelVariation = clone_JSON_O(currentModel);
						modifyBindingMode(modelVariation,  deltaLeft,  deltaRight, deltaFlank, dIntMax, componentList);
						newVariations.put(testLabel, modelVariation);
					}
				} 
				
				if(newVariations.size()>0)
					return newVariations;
			} 
			
			//2.2 Determines if the optimizer should move on to the next fitting stage.
			if(fittingStages==null)
				break;
			else {
				if ( (iFittingStage+1)<fittingStages.length() ) {
					setFittingStage(iFittingStage+1);
					continue;
				} else 
					break;
			}
		} 

		System.out.println("> No varitions possible for "+componentName+".");
		//No variation possible.
		variationsOptimized = true;
		return null;
		
	}
	
	private double computeBaseInformation(int x) {
		
		//Computes the probability of bases
		double[] probability = new double[nMono];
		double expSum        = 0;
		for(int i=0; i<nMono; i++) {
			double temp      = Math.exp(monoBetas[nMono*x+i]);
			expSum          += temp;
			probability[i]   = temp;
		}
		for(int i=0; i<nMono; i++)
			probability[i]  /= expSum;
		
		//Computes the information
		double log2          = Math.log(2);
		double information   = Math.log(letterOrder.length())/log2;		
		for(int i=0; i<nMono; i++) 
			if(probability[i]>0)
				information     += probability[i] * Math.log(probability[i]) / log2;
		
		return information; 
		
	}
	
	
	@Override
	public void indicateModelChoice(String bestKey) {
		if(bestKey!=null) {
			currentShift  = Integer.parseInt(bestKey.split(";")[2]);
			variationName = bestKey;
		}
		return;
	}
	
	@Override
	public void seed_component(JSONObject config) {
		
		double[] monoSeed = null;
		
		String coefficientKey = "modelSeeding";
		
		if(config.has(coefficientKey) && 
				config.getJSONObject(coefficientKey).has(componentKey) && 
				config.getJSONObject(coefficientKey).getJSONArray(componentKey).length() > iComp) {
			
			JSONObject oBm = config.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
			
			//Sets the seed scale
			double seedScale = oBm.has("seedScale") ? oBm.getDouble("seedScale") : 1.0 ;
			
			if(oBm.has("mononucleotideIUPAC")) {
				
				if(letterOrder.length()>=4 && letterOrder.subSequence(0, 4).equals("ACGT")) {
					String iupacString = oBm.getString("mononucleotideIUPAC");

					if(iupacString.length() != k)
						throw new IllegalArgumentException("The length of the IUPAC "+iupacString+" does not match k="+k);
					monoSeed = new double[nMono*k];
					
					for(int i=0; i<iupacString.length(); i++) {
						double[] r = null;
						double m=seedScale;
						switch(iupacString.toUpperCase().charAt(i)) {
							case 'A': r = new double[]{ 0,-m,-m,-m};break;
							case 'C': r = new double[]{-m, 0,-m,-m};break;
							case 'G': r = new double[]{-m,-m, 0,-m};break;
							case 'T': r = new double[]{-m,-m,-m, 0};break;
							case 'R': r = new double[]{ 0,-m, 0,-m};break;
							case 'Y': r = new double[]{-m, 0,-m, 0};break;
							case 'S': r = new double[]{-m, 0, 0,-m};break;
							case 'W': r = new double[]{ 0,-m,-m, 0};break;
							case 'K': r = new double[]{-m,-m, 0, 0};break;
							case 'M': r = new double[]{ 0, 0,-m,-m};break;
							case 'B': r = new double[]{-m, 0, 0, 0};break;
							case 'D': r = new double[]{ 0,-m, 0, 0};break;
							case 'H': r = new double[]{ 0, 0,-m, 0};break;
							case 'V': r = new double[]{ 0, 0, 0,-m};break;
							case 'N': r = new double[]{ 0, 0, 0, 0};break;
							default: throw new IllegalArgumentException("Invalid IUPAC code: "+iupacString.charAt(i));
						}
						
						for(int n=0; n<nMono; n++)
							monoSeed[nMono*i+n] = n<4 ? r[n] : -m;
							
					}
				} else {
					throw new IllegalArgumentException("IUPAC code seed can only be used when the first four letters in the ordered alphabet are 'ACGT'");
				}
			} else if(oBm.has("mononucleotideString")) {
				String monoString = oBm.getString("mononucleotideString");
				
				if(monoString.length() != k)
					throw new IllegalArgumentException("The length of the seeding string "+monoString+" does not match k="+k);
				monoSeed = new double[nMono*k];
				
				for(int i=0; i<monoString.length(); i++) {
					Character c = monoString.charAt(i);
					double m    = seedScale;
					int iLet    = letterOrder.indexOf(c);
					if(iLet>=0) {
						//Puts all non-matching letters to -m.
						for(int iM=0; iM<nMono; iM++) {
							if(iM!=iLet) 
								monoSeed[nMono*i+iM] = -m;
						}
						//monoSeed[nMono*i+iLet] =  m;
					} else if(c=='.') {
						//Do nothing - all should be zero
					} else {
						throw new IllegalArgumentException("Seeding letter '"+c+"' not in alphabet '"+letterOrder+"'.");
					}
				}
			}
		}
		
		if(monoSeed==null && monoBetas!=null)
			monoSeed = random_d(monoBetas, 0.01);
		
		monoBetas = monoSeed;
	}

	@Override
	int packModel_component(JSONObject packing, int iFirst) {
		
		String coefficientKey = "packing";
		addEmptyJSON_component_O(packing, coefficientKey, componentKey, iComp);

		int iCurr = iFirst;
		Integer iLast;
		
		//Builds packing indices for the binding mode symmetry
		
		ArrayList<int[]> symmetryStringPacking = null;
		if( (!symmetryString.equals("null")) && k>0) {
				BindingModeSymmetry symmetry = new BindingModeSymmetry(symmetryString, k, letterOrder, letterComplement);
				symmetryStringPacking        = symmetry.constructPackingVectors(iCurr+1);
		} 
		
		if(fitMono) {
			
			iLast                                           = null;
			int[] monoRange                                 = null;
			ArrayList<int[]> diRange                        = null;
			ArrayList<int[]> monoModificationRange          = null;
			ArrayList<ArrayList<int[]>> diModificationRange = null;

			//Packs mononucleotides
			///////////////////////
			
			//Mononucleotide betas
			if(symmetryStringPacking==null) {
				monoRange = range_d(monoBetas,          iCurr+1);
				iLast = last_i(monoRange);
			} else {
				monoRange = symmetryStringPacking.get(0);
				if(monoRange!=null)
					for(int j=0; j<monoRange.length; j++)
						if( iLast==null || monoRange[j] > iLast)
							iLast = monoRange[j];
			}
			iCurr = iLast==null ? iCurr : iLast;
			saveToJSON_mono_i(packing, coefficientKey, monoRange);
			
			//Mononucleotide beta modifications
			if(modifications.size()>0) {
				//Naive modifications.
				monoModificationRange = range_Ad(monoBetaShiftValues, iCurr+1);
				iLast = last_Ai(monoModificationRange);
				iCurr = iLast==null ? iCurr : iLast;
				
				//Accounts for symmetry
				for(int iMod=0; iMod<monoBetaShiftIndices.size(); iMod++){
					HashMap<Integer,Integer> featureRepresentative = new HashMap<Integer,Integer>();
					int[] monoModI                = monoBetaShiftIndices.get(iMod);
					int[] monoModPackedI          = monoModificationRange.get(iMod);

					for(int iFeat=0; iFeat<monoModPackedI.length; iFeat++) {
						int naiveIndex           = monoModPackedI[iFeat];
						int packedMonoIndex      = monoRange[monoModI[iFeat]];
						if(featureRepresentative.containsKey(packedMonoIndex))
							monoModPackedI[iFeat] = featureRepresentative.get(packedMonoIndex);
						else
							featureRepresentative.put(packedMonoIndex, naiveIndex);
					}
				}
			}

			if(fitDi) {

				//Dinucleotide betas
				iLast = null;
				if(symmetryStringPacking==null) {
					diRange= range_Ad(diBetas,            iCurr+1);
					iLast = last_Ai(diRange);
				} else {
					diRange = new ArrayList<int[]>();
					for(int iD=0; iD<diBetas.size(); iD++) {
						int[] currentPacking = symmetryStringPacking.get(iD+1);
						diRange.add(currentPacking);
						//Finds the largest index
						if(currentPacking!=null)
							for(int j=0; j<currentPacking.length; j++)
								if( iLast==null || currentPacking[j] > iLast)
									iLast = currentPacking[j];
					}
				}
				iCurr = iLast==null ? iCurr : iLast;
				saveToJSON_di_Ai(packing, coefficientKey, diRange);
				
				//Dinucleotide beta modifications
				if(modifications.size()>0) {
					
					//Naive modifications
					diModificationRange = range_AAd(diBetaShiftValues, iCurr+1);
					iLast = last_AAi(diModificationRange);
					iCurr = iLast==null ? iCurr : iLast;

					//Accounts for symmetry
					for(int iMod=0; iMod<diModificationRange.size(); iMod++){
						HashMap<Integer,Integer> featureRepresentative = new HashMap<Integer,Integer>();
						for(int iD=0; iD<diModificationRange.get(iMod).size(); iD++) {
							int[] diModI             = diBetaShiftIndices.get(iMod).get(iD);
							int[] diModPackedI       = diModificationRange.get(iMod).get(iD);

							for(int iFeat=0; iFeat<diModI.length; iFeat++) {

								int naiveIndex          = diModPackedI[iFeat];
								int packedDiIndex       = diRange.get(iD)[diModI[iFeat]];
								if(featureRepresentative.containsKey(packedDiIndex))
									diModPackedI[iFeat] = featureRepresentative.get(packedDiIndex);
								else
									featureRepresentative.put(packedDiIndex, naiveIndex);
							}
						}
					}
				}
			}
			
			//Saves the modification indices. 
			if(monoModificationRange!=null || diModificationRange!=null)
				saveToJSON_modification_Ai_AAi(packing, coefficientKey, monoModificationRange, diModificationRange);
		}
		

		if(fitPositionBias) {
			
			//Creates new index list assuming no constraints
			ArrayList<ArrayList<int[]>> positionBiasRange = range_AAd(positionBiasBetas, iCurr+1);
			iLast = last_AAi(positionBiasRange);
			iCurr = iLast==null ? iCurr : iLast;
			
			//Imposes singleStrand
			if(singleStrand)
				for(int iExp=0; iExp<positionBiasRange.size(); iExp++) {
					int[] v = positionBiasRange.get(iExp).get(1);
					for(int iX=0; iX<v.length; iX++)
						v[iX] = -1;
				}
			
			//Imposes position-bias bins:
			if(wPositionBiasBin!=1)
				for(int iExp=0; iExp<positionBiasRange.size(); iExp++) {
					for(int iS=0; iS<2; iS++) {
						int[] v = positionBiasRange.get(iExp).get(iS);
						for(int iX=0; iX<v.length; iX++)
							v[iX] = v[iX - (iX%wPositionBiasBin)];
					}
				}
			
			//Imposes experimentSpecificPositionBias=False
			if(!experimentSpecificPositionBias)
				for(int iExp=0; iExp<positionBiasRange.size(); iExp++) {
					for(int iS=0; iS<2; iS++) {
						int[] v0 = positionBiasRange.get(0).get(iS);
						int[] v1 = positionBiasRange.get(iExp).get(iS);
						if(v0.length != v1.length)
				        	throw new java.lang.RuntimeException("ERROR: All experiments must have the same number of binding frames if BindingMode.experimentSpecificPositionBias=False.");

						for(int iX=0; iX<v1.length; iX++)
							v1[iX] = v0[iX];

					}
				}
			
			saveToJSON_positionBias_AAi(packing, coefficientKey, positionBiasRange);
			
		}

		if(fitActivity) {
			
			//Setting up activities
			ArrayList<int[]> activityRange = range_Ad(activityBetas, iCurr+1);
			iLast        = last_Ai(activityRange);
			iCurr        = iLast==null ? iCurr : iLast;
			int nExps    = activityRange.size();
			int nRounds0 = activityRange.get(0).length;
			
			//Removes round-specific activities if appropriate. 
			if(roundSpecificActivity==false) 
				for(int iExp=0; iExp<nExps; iExp++)
					for(int iRound=0; iRound<activityRange.get(iExp).length; iRound++)
						activityRange.get(iExp)[iRound] = activityRange.get(iExp)[0]; 
			//Removes experiment-specific activities if appropriate.
			if(experimentSpecificActivity==false && nExps>1) {
				for(int iExp=1; iExp<nExps; iExp++) {
					if(activityRange.get(iExp).length != nRounds0)
			        	throw new java.lang.RuntimeException("ERROR: All experiments must have the rounds if BindingMode.experimentSpecificActivity=False.");

					for(int iRound=0; iRound<nRounds0; iRound++) {
						activityRange.get(iExp)[iRound] = activityRange.get(0)[iRound];
					}
				}
			}

			saveToJSON_activity_Ai(packing, coefficientKey, activityRange);
			
		}
		
		return iCurr;

	}

	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {

		addEmptyJSON_component_O(in,        coefficientKey, componentKey, iComp);
		saveToJSON_mono_d(in,               coefficientKey, zero_d(monoBetas));
		saveToJSON_di_Ad(in,                coefficientKey, zero_Ad(diBetas));
		if(positionBiasBetas!=null) 
			saveToJSON_positionBias_AAd(in, coefficientKey, zero_AAd(positionBiasBetas));
		saveToJSON_activity_Ad(in,          coefficientKey, zero_Ad(activityBetas));
		
		saveToJSON_modification_Ad_AAd (in, coefficientKey, zero_Ad(monoBetaShiftValues), zero_AAd(diBetaShiftValues));
		
		return;
	}

	@Override
	void saveToJSON_settings(JSONObject out) {

		String coefficientKey = "modelSettings";
		addEmptyJSON_component_O(out,       coefficientKey, componentKey, iComp);
		JSONObject oBm   =                  out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBm.put("size",                     k);
		oBm.put("flankLength",              flankLength);
		oBm.put("dinucleotideDistance",     dIntMax);
		oBm.put("positionBias",             usePositionBias);
		oBm.put("singleStrand",             singleStrand);
		oBm.put("fitLogActivity",           fitLogActivity);
		
		JSONArray aMod = new JSONArray();
		for(int iMod=0; iMod<modifications.size(); iMod++) {
			JSONObject oMod = new JSONObject();
			oMod.put("name", modifications.get(iMod));
			
			//Saves mononucleotides.
			JSONArray aMonoMods = new JSONArray();
			Character[] mc = monoMod.get(iMod);
			for(int iMono=0; iMono < mc.length; iMono++)
				aMonoMods.put(mc[iMono].toString());
			oMod.put("mononucleotide", aMonoMods);

			//Saves dinucleotides
			JSONArray aDiMods = new JSONArray();
			for(int d=0; d<dIntMax; d++) {
				JSONArray aDiModsI = new JSONArray();
				diMod.get(iMod);
				diMod.get(iMod).get(d);
				Character[][] dc = diMod.get(iMod).get(d);
				for(int iDi=0; iDi<dc.length; iDi++)
					aDiModsI.put(""+dc[iDi][0]+dc[iDi][1]);
				aDiMods.put(aDiModsI);
			}
			oMod.put("dinucleotide", aDiMods);
			aMod.put(oMod);
		}
		oBm.put("modifications", aMod);
		
	}

	@Override
	void readFromJSON_settings(JSONObject in) {
		
		String coefficientKey = "modelSettings";
		
		JSONObject oBm = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		k              = oBm.getInt("size");
		flankLength    = oBm.getInt("flankLength");
		dIntMax        = oBm.getInt("dinucleotideDistance");
		usePositionBias= oBm.getBoolean("positionBias");
		singleStrand   = oBm.getBoolean("singleStrand");
		fitLogActivity = oBm.getBoolean("fitLogActivity");
		dInt           = Math.max(0,Math.min(k-1, dIntMax));
		
		//Loads modifications
		JSONArray aMod = oBm.getJSONArray("modifications");
		modifications        = new ArrayList<String> ();
		monoMod              = new ArrayList<Character[]>();
		diMod                = new ArrayList<ArrayList<Character[][]>>();
		monoBetaShiftIndices = new ArrayList<int[]>();
		diBetaShiftIndices   = new ArrayList<ArrayList<int[]>>();

		for(int iMod=0; iMod<aMod.length(); iMod++) {
			JSONObject oMod = aMod.getJSONObject(iMod);
			String modName  = oMod.getString("name");
			modifications.add(modName);
			
			//Loading mononucleotide modifications and building indices
			JSONArray aMono                = oMod.getJSONArray("mononucleotide");
			int nMonoMod                   = aMono.length();
			Character[] monoChars          = new Character[nMonoMod];
			int[] newMonoIndices           = new int[k*nMonoMod];
			for(int iMono=0; iMono<nMonoMod; iMono++) {
				monoChars[iMono] = aMono.getString(iMono).charAt(0);
				for(int x=0; x<k; x++)
							newMonoIndices[x*nMonoMod+iMono] = nMono*x + letterOrder.indexOf(monoChars[iMono]);
			}
			monoMod.add(monoChars);
			monoBetaShiftIndices.add(newMonoIndices);

			//Loading dinucleotide
			JSONArray aDi                     = oMod.getJSONArray("dinucleotide");
			ArrayList<Character[][]> newDiMod = new ArrayList<Character[][]>();
			ArrayList<int[]> newDI_Ai         = new ArrayList<int[]>(); 
			for(int d=0; d<dIntMax; d++) { //Loops over spacing
				//Creating a list of modifications
				Character[][] newDiChars;
				if(d<aDi.length()) {
					JSONArray aDiI = aDi.getJSONArray(d);
					newDiChars = new Character[aDiI.length()][];
					for(int j=0; j<aDiI.length(); j++) {
						newDiChars[j] = new Character[2];
						newDiChars[j][0] = aDiI.getString(j).charAt(0);
						newDiChars[j][1] = aDiI.getString(j).charAt(1);
					}
				} else {
					newDiChars = new Character[0][0];
				}
				newDiMod.add(newDiChars);
				
				//Creating the list of indices
				int nDiMod       = newDiChars.length;
				int[] newDI_i = new int[(k-d-1)*nDiMod];
				for(int iDi=0; iDi<nDiMod; iDi++) { //Loops over dinucleotides.
					int diShifts = nMono * letterOrder.indexOf(newDiChars[iDi][0]) +
							letterOrder.indexOf(newDiChars[iDi][1]);
					for(int x=0; x<k-d-1; x++) //Loops over positions.
						newDI_i[x*nDiMod+iDi] = nDi*x + diShifts;
				}
				newDI_Ai.add(newDI_i);
			}		
			diMod.add(newDiMod);
			diBetaShiftIndices.add(newDI_Ai);

		}
		
		allocateParameters();

	}
	
	@Override
	void saveToJSON_constraints(JSONObject out) {
		
		String coefficientKey = "modelFittingConstraints";
		addEmptyJSON_component_O(out, coefficientKey, componentKey, iComp);
		
		JSONObject oBm        =                   out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);

		oBm.put("positionBiasBinWidth",           wPositionBiasBin);
		oBm.put("optimizeSize",                   fitK);
		oBm.put("optimizeSizeHeuristic",          optimizeSizeHeuristic);
		oBm.put("optimizeFlankLength",            fitFlankLength);
		oBm.put("optimizeMotifShift",             shiftMotif);
		oBm.put("optimizeMotifShiftHeuristic",    heuristicShift);
		if(fittingStages!=null)
			oBm.put("fittingStages",              ModelComponent.clone_JSON_A(fittingStages));
		oBm.put("minSize",                        minSize);
		oBm.put("maxSize",                        maxSize);
		oBm.put("maxFlankLength",                 maxFlankLength);
		oBm.put("symmetryString",                 symmetryString);
		oBm.put("roundSpecificActivity",          roundSpecificActivity);
		oBm.put("experimentSpecificPositionBias", experimentSpecificPositionBias);
		oBm.put("experimentSpecificActivity",     experimentSpecificActivity);
		oBm.put("informationThreshold",           informationThreshold);
		
		
	}

	@Override
	void readFromJSON_constraints(JSONObject in) {
		if(in.has("modelFittingConstraints")) {
			JSONObject oOpt       = in.getJSONObject("modelFittingConstraints").getJSONArray(componentKey).getJSONObject(iComp);
			wPositionBiasBin      = oOpt.getInt(     "positionBiasBinWidth");
			
			fitK                  = oOpt.getBoolean( "optimizeSize");
			fitFlankLength        = oOpt.getBoolean( "optimizeFlankLength");
			optimizeSizeHeuristic = oOpt.has("optimizeSizeHeuristic")       ? oOpt.getBoolean( "optimizeSizeHeuristic") : false;
			shiftMotif            = oOpt.getBoolean( "optimizeMotifShift");
			heuristicShift        = oOpt.has("optimizeMotifShiftHeuristic") ? oOpt.getBoolean( "optimizeMotifShiftHeuristic") : false;
			if(oOpt.has("fittingStages")) {
				fittingStages     = ModelComponent.clone_JSON_A(oOpt.getJSONArray("fittingStages"));
				if(fittingStages.length()>0)
					setFittingStage(iFittingStage);
			}
			
			minSize               = oOpt.has("minSize")                     ? oOpt.getInt("minSize")        : -1;
			maxSize               = oOpt.has("maxSize")                     ? oOpt.getInt("maxSize")        : -1;
			maxFlankLength        = oOpt.has("maxFlankLength")              ? oOpt.getInt("maxFlankLength") : -1;

			symmetryString        = oOpt.getString(  "symmetryString");
			roundSpecificActivity = oOpt.getBoolean( "roundSpecificActivity");
			experimentSpecificPositionBias
			                      = oOpt.getBoolean( "experimentSpecificPositionBias");
			experimentSpecificActivity
			                      = oOpt.getBoolean( "experimentSpecificActivity");
			if(oOpt.has("informationThreshold"))
				informationThreshold = oOpt.getDouble( "informationThreshold");
			else
				informationThreshold = 0.1; //Note: This could be removed (Only here to be back-compatible)
			
		}
	}

	@Override
	public void regularizationToJSON_component(JSONObject out, String coefficientKey) {

		addEmptyJSON_component_O(out,        coefficientKey, componentKey, iComp);

		//Penalization of PSAM
		if(fitMono && k>0) {
			double[] monoReg = zero_d(monoBetas);
			double meanValue;
			for(int x=0; x<k; x++) {
				// Computes mean value of mono.
				meanValue=0;
				for(int n=0; n<nMono; n++) 
					meanValue 													  += monoBetas[nMono*x+n]/nMono;
				// Adds mononucleotide error contribution.
				for(int n=0; n<nMono; n++)
					monoReg[nMono*x+n] 								               = monoBetas[nMono*x+n] - meanValue; 
			}
			saveToJSON_mono_d(out, coefficientKey, monoReg);			
		}
		
		//Penalization of dinucleotide PSAM
		if(fitDi && diBetas != null ) {
			ArrayList<double[]> diReg = zero_Ad(diBetas);
			double meanValue;
			for(int d=1; d<=dInt; d++){
				for(int x=0; x<k-d; x++) {
					// Computes the mean value for dinucleotide
					meanValue=0;
					for(int n=0; n<nDi; n++) 
						meanValue												  += diBetas.get(d-1)[nDi*x+n] / nDi;
					// Adds mononucleotide error contribution.
					for(int n=0; n<nDi; n++) 
						diReg.get(d-1)[nDi*x+n]                                    = diBetas.get(d-1)[nDi*x+n] - meanValue; 
				}
			}
			saveToJSON_di_Ad(out, coefficientKey, diReg);			
		}


		//Penalization of position bias
		if( usePositionBias && k>0 ){
			ArrayList<ArrayList<double[]>> positionBiasReg = zero_AAd(positionBiasBetas);
			for(int iExp=0; iExp<positionBiasBetas.size(); iExp++)
				for(int s=0; s<2; s++)
					for(int x=0; x<positionBiasBetas.get(iExp).get(s).length; x++)
						positionBiasReg.get(iExp).get(s)[x]                        = positionBiasBetas.get(iExp).get(s)[x];
			saveToJSON_positionBias_AAd(out, coefficientKey, positionBiasReg);
		}
		
		if(fitLogActivity)
			saveToJSON_activity_Ad(out,          coefficientKey, zero_Ad(activityBetas));
		else {
			saveToJSON_activity_Ad(out,          coefficientKey, activityAlphas);
		}
		
		saveToJSON_modification_Ad_AAd (out,     coefficientKey, monoBetaShiftValues, diBetaShiftValues);
				
		return;
		
	}

	@Override
	void saveToJSON_parameters(JSONObject out, String coefficientKey) {

		addEmptyJSON_component_O(out,        coefficientKey, componentKey, iComp);
		saveToJSON_mono_d(out,               coefficientKey, monoBetas);
		saveToJSON_di_Ad(out,                coefficientKey, diBetas);

		if(positionBiasBetas!=null) 
			saveToJSON_positionBias_AAd(out, coefficientKey, positionBiasBetas);

		if(fitLogActivity)
			saveToJSON_activity_Ad(out,          coefficientKey, activityBetas);
		else {
			saveToJSON_activity_Ad(out,          coefficientKey, activityAlphas);
		}
		
		saveToJSON_modification_Ad_AAd(out,  coefficientKey, monoBetaShiftValues, diBetaShiftValues);
		
		return;
		
	}

	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {

		JSONObject oBM        = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		monoBetas             = readFromJSON_d(oBM.getJSONArray( "mononucleotide"));
		diBetas               = readFromJSON_Ad(oBM.getJSONArray("dinucleotide"));
		if(usePositionBias && k>0)
			positionBiasBetas = readFromJSON_AAd(oBM.getJSONArray("positionBias"));
		if(fitLogActivity)
			activityBetas     = readFromJSON_Ad(oBM.getJSONArray("activity"));
		else
			activityAlphas    = readFromJSON_Ad(oBM.getJSONArray("activity"));
		

		//Reads modifications.
		monoBetaShiftValues   = new ArrayList<double[]>();
		diBetaShiftValues     = new ArrayList<ArrayList<double[]>>();
		JSONArray aMod        = oBM.getJSONArray("modifications");
		for(int iMod=0; iMod<aMod.length(); iMod++ ) {
			JSONObject oMod = aMod.getJSONObject(iMod);
			monoBetaShiftValues.add(readFromJSON_d(oMod.getJSONArray("mononucleotide")));
			diBetaShiftValues.add(readFromJSON_Ad(oMod.getJSONArray("dinucleotide")));
		}
				
		updateAlphas();
		//updateSlidingWindow();
	}

	//Methods for saving double objects
	void saveToJSON_mono_d(JSONObject in, String coefficientKey, double[] monoBetas) {
		JSONObject oBM = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("mononucleotide", JSONArrayConvert_d(monoBetas));
		return;
	}
	
	void saveToJSON_di_Ad(JSONObject in, String coefficientKey, ArrayList<double[]> diBetas) {
		JSONObject oBM = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("dinucleotide", JSONArrayConvert_Ad(diBetas));
		return;
	}
	
	void saveToJSON_positionBias_AAd(JSONObject in, String coefficientKey, ArrayList<ArrayList<double[]>> positionBiasBetas) {
		JSONObject oBM = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("positionBias", JSONArrayConvert_AAd(positionBiasBetas));
		return;
	}
	
	void saveToJSON_positionBias_Ad(JSONObject in, String coefficientKey, ArrayList<double[]> positionBiasBetas, int iExp) {
		JSONObject oBM = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);

		if(positionBiasBetas != null ){
			
			//Gets the positionBias object.
			JSONArray aPB;
			if(!oBM.has("positionBias"))
				oBM.put("positionBias", new JSONArray());
			aPB = oBM.getJSONArray("positionBias");
			
			//Makes sure all he
			while(aPB.length() < countTables.size())
				aPB.put(new JSONArray());
			
			//Saves the specific betas.
			aPB.put(iExp, JSONArrayConvert_Ad(positionBiasBetas));
		
		}


		return;
	}
	
	void saveToJSON_activity_Ad(JSONObject out, String coefficientKey, ArrayList<double[]> activityBetas) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("activity", JSONArrayConvert_Ad(activityBetas));
			
		return;
	}
	
	void saveToJSON_activity_d(JSONObject out, String coefficientKey, double[] activityBetas, int iExp) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		
		if(!oBM.has("activity"))
			oBM.put("activity", new JSONArray());
		JSONArray aA = oBM.getJSONArray("activity");
		
		//Makes sure all he
		while(aA.length() < countTables.size())
			aA.put(new JSONArray());
		
		//Saves the specific betas.
		aA.put(iExp, JSONArrayConvert_d(activityBetas));

		return;
	}
	
	void saveToJSON_modification_Ad_AAd(JSONObject out, String coefficientKey, 
			ArrayList<double[]> monoModificationShifts, ArrayList<ArrayList<double[]>> diModificationShifts) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);

		if(monoModificationShifts!=null || diModificationShifts!=null) {

			if(!oBM.has("modifications"))
				oBM.put("modifications", new JSONArray());
			JSONArray aMod = oBM.getJSONArray("modifications");
			
			//Making sure all objects exist.
			while(aMod.length() < modifications.size())
				aMod.put(new JSONObject());

			//Saves the mononucleotide and dinucleotide data
			for(int iMod=0; iMod<modifications.size(); iMod++) {
				JSONObject oMod = aMod.getJSONObject(iMod);
				if(monoModificationShifts!=null)
					oMod.put("mononucleotide", JSONArrayConvert_d(monoModificationShifts.get(iMod)));
				if(diModificationShifts != null)
					oMod.put("dinucleotide",   JSONArrayConvert_Ad(diModificationShifts.get(iMod)));
				
				//aMod.put(iMod, oMod);
			}
		}
	}
	
	//Methods for saving long objects
	void saveToJSON_mono_i(JSONObject in, String coefficientKey, int[] monoBetas) {
		JSONObject oBM = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("mononucleotide", JSONArrayConvert_i(monoBetas));
		return;
	}
	
	void saveToJSON_di_Ai(JSONObject in, String coefficientKey, ArrayList<int[]> diBetas) {
		JSONObject oBM = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("dinucleotide", JSONArrayConvert_Ai(diBetas));
		return;
	}
	
	void saveToJSON_positionBias_AAi(JSONObject in, String coefficientKey, ArrayList<ArrayList<int[]>> positionBiasBetas) {
		JSONObject oBM = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("positionBias", JSONArrayConvert_AAi(positionBiasBetas));
		return;
	}
	
	void saveToJSON_activity_Ai(JSONObject out, String coefficientKey, ArrayList<int[]> activityBetas) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("activity", JSONArrayConvert_Ai(activityBetas));
		return;
	}
	
	void saveToJSON_modification_Ai_AAi(JSONObject out, String coefficientKey, 
			ArrayList<int[]> monoModificationRange, ArrayList<ArrayList<int[]>> diModificationRange) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		
		if(!oBM.has("modifications"))
			oBM.put("modifications", new JSONArray());
		JSONArray aMod = oBM.getJSONArray("modifications");
		
		//Making sure all objects exist.
		while(aMod.length() < modifications.size())
			aMod.put(new JSONObject());

		//Saves the mononucleotide and dinucleotide data
		for(int iMod=0; iMod<modifications.size(); iMod++) {
			JSONObject oMod = aMod.getJSONObject(iMod);
			if(monoModificationRange!=null)
				oMod.put("mononucleotide", JSONArrayConvert_i(monoModificationRange.get(iMod)));
			if(diModificationRange!=null)
				oMod.put("dinucleotide", JSONArrayConvert_Ai(diModificationRange.get(iMod)));
			//aMod.put(iMod, oMod);
		}
	}
	
	public double computeInformationWeightedMean(){
		//Running sums used to compute the weighted mean position in the binding model.
		double wSum = 0, wxSum = 0;
				
		//Loops across the PSAM.
		for(int i=0; i<k; i++) {
			
			//Computes the probability of each base.
			double[] z = new double[nMono], p = new double[nMono];
			for(int c=0; c<nMono; c++)
				z[c] = Math.exp(monoBetas[nMono*i+c]);
			double zSum = Array.sum(z);
			for(int c=0; c<nMono; c++)
				p[c] = z[c] / zSum;

			double information = 0;
			for(int c=0; c<nMono; c++)
				information += p[c] * Math.log(p[c]*nMono) / Math.log(2);
			
			wSum  += information;
			wxSum += information * i;
		}
		
		return wxSum / wSum;
		
	}
	
	
	public double computeVarianceWeightedMean(){

		//Running sums uset to compute the weighted mean position in the binding model.
		double wSum = 0, wxSum = 0;
		
		//Conditional mean energy.
		double[] expectedE = new double[nMono];
		
		//Loops across the PSAM.
		for(int i=0; i<k; i++) {
			//Resets the deviation vector 
			for(int c=0; c<nMono; c++) expectedE[c] = 0;
			
			//Computes the computational mean energy.
			for(int c=0; c<nMono; c++) {
				for(int j=0; j<k; j++){
					for(int d=0; d<nMono; d++){
						if( j<i && i-j<=dIntMax )
							//(j,d) before (i,c)
							expectedE[c] += diBetas.get(i-j-1)[ nDi*j + nMono*d + c ] / nMono;
							//(j,d) after (i,c)
						else if( j>i && j-i<=dIntMax )
							expectedE[c] += diBetas.get(j-i-1)[ nDi*i + nMono*c + d ] / nMono;
						else if( i==j )
							expectedE[c] += monoBetas[nMono*i+c];
					}
				}
			}
			
			//Computes the deviation from the mean conditional mean and then the variance
			double meanExpectedE = Array.mean(expectedE), var=0;
			for(int c=0; c<nMono; c++) {
				var += Math.pow(expectedE[c]-meanExpectedE, 2) / nMono;
			}	
			wSum  += var;
			wxSum += var * i;
		}
		
		return wxSum / wSum;
	}
	
	
	static public void printJSONObjectCoefficients(JSONObject model, String coefficientKey, int iBM) {
		
		JSONObject oBM = model.getJSONObject(coefficientKey).getJSONArray("bindingModes").getJSONObject(iBM);
		
		
		System.out.println("Mononucleotide:    "
				+Misc.formatVector_d(readFromJSON_d(oBM.getJSONArray("mononucleotide"))));
		if(oBM.has("dinucleotide")) {
			JSONArray aDI = oBM.getJSONArray("dinucleotide");
			for(int iD=0; iD<aDI.length(); iD++) {
				System.out.println("Dinucleotide(d="+(iD+1)+"): "
						+Misc.formatVector_d(readFromJSON_d(aDI.getJSONArray(iD))));
			}
			
		}

		if(oBM.has("activity")) {
			JSONArray aAct = oBM.getJSONArray("activity");
			for(int iExp=0; iExp<aAct.length(); iExp++) {
				System.out.println("Activity(exp="+iExp+"):   "
						+Misc.formatVector_d(readFromJSON_d(aAct.getJSONArray(iExp))));
			}
		}

		if(oBM.has("positionBias")) {
			JSONArray aPB = oBM.getJSONArray("positionBias");
			for(int iExp=0; iExp<aPB.length(); iExp++) {
				for(int s=0; s<2; s++) {
					System.out.println("Position Bias(exp="+iExp+", s="+s+"): "
							+Misc.formatVector_d(readFromJSON_d(aPB.getJSONArray(iExp).getJSONArray(s))));
				}
			}
		}
		
		if(oBM.has("modifications")) {
			JSONArray aMod = oBM.getJSONArray("modifications");
			for(int iMod=0; iMod<aMod.length(); iMod++) {
				JSONObject oMod = aMod.getJSONObject(iMod);
				System.out.println("Modification "+iMod+".");
				JSONArray aMono = oMod.getJSONArray("mononucleotide");
				if(aMono.length()>0)
					System.out.println("Mononucleotide:"+Misc.formatVector_d(readFromJSON_d(aMono)));
				JSONArray aDi = oMod.getJSONArray("dinucleotide");
				if(aDi.length()>0) {
					System.out.print("Dinucleotide ("+iMod+"):");
					for(int d=0; d<aDi.length(); d++)
						System.out.print(Misc.formatVector_d(readFromJSON_d(aDi.getJSONArray(d)))+",");
					System.out.println("");
				}

			}
		}
	}
	
	
	public void modifyBindingMode(JSONObject model, int deltaLeft, int deltaRight, int deltaFlankLength, int newDIntMax, ArrayList<ModelComponent> componentList) {
		
		//Identifies the JSON output object. 
		JSONObject oBMSetting  = model.getJSONObject("modelSettings").getJSONArray("bindingModes").getJSONObject(iComp);
		JSONObject oBMCoeff    = model.getJSONObject("coefficients" ).getJSONArray("bindingModes").getJSONObject(iComp);

		//Updates binding mode settings 
		int newK           = k + deltaLeft + deltaRight;
		int newFlankLength = flankLength + deltaFlankLength;
		int oldDInt        = dInt;
		int newDInt        = Math.max(0,Math.min(newK-1, newDIntMax));
		oBMSetting.put("size", newK);
		oBMSetting.put("flankLength", newFlankLength);
		oBMSetting.put("dinucleotideDistance", newDIntMax);
		
		//Updates nucleotide betas
		oBMCoeff.put("mononucleotide", JSONArrayConvert_d(Array.padArray(monoBetas, nMono*deltaLeft, nMono*deltaRight)));
		
		//Updates dinucleotide betas
		ArrayList<double[]> newDiBetas = new ArrayList<double[]>();  
		for (int d=1; d <= newDInt; d++)
			newDiBetas.add(d<=oldDInt ? Array.padArray(diBetas.get(d-1), nDi*deltaLeft, nDi*deltaRight) : new double[nDi*(newK-d)]);
		oBMCoeff.put("dinucleotide", JSONArrayConvert_Ad(newDiBetas));
		
		//Updates the position bias vector
		if(usePositionBias){
			if(newK>0) {

				for(int iExp=0; iExp<positionBiasBetas.size(); iExp++ ){
					for(int s1=0; s1<2; s1++) {
						oBMCoeff.getJSONArray("positionBias").getJSONArray(iExp).put(s1, 
							JSONArrayConvert_d( k>0 
								? Array.padArray(positionBiasBetas.get(iExp).get(0), -deltaLeft  + deltaFlankLength, -deltaRight + deltaFlankLength)
								: new double[countTables.get(iExp).l + 2*newFlankLength - newK + 1]));
					}
				}
			}
		} 
		
		//Shifting binding mode modifications
		for(int iMod=0; iMod<modifications.size(); iMod++) {
			JSONObject oMod = oBMCoeff.getJSONArray("modifications").getJSONObject(iMod);
			//Updates mononucleotide betas
			int nMonoMod = monoMod.get(iMod).length;
			oMod.put("mononucleotide", JSONArrayConvert_d(Array.padArray(monoBetaShiftValues.get(iMod), nMonoMod*deltaLeft, nMonoMod*deltaRight)));

			//Updates dinucleotide betas
			ArrayList<double[]> newDiModBetas_Ad = new ArrayList<double[]>(); 
			for (int d=1; d <= newDInt; d++) {
				int nDiMod          = diMod.get(iMod).get(d-1).length;
				newDiModBetas_Ad.add( d<=oldDInt ? Array.padArray(diBetaShiftValues.get(iMod).get(d-1), nDiMod*deltaLeft,   nDiMod*deltaRight) : new double[nDiMod*(newK-d)] );
			}
			oMod.put("dinucleotide",   JSONArrayConvert_Ad(newDiModBetas_Ad));
		}
		
		//Updates the interactions.
		for( ModelComponent comp : componentList )
			if ( comp instanceof BindingModeInteraction )
				((BindingModeInteraction) comp).modifyBindingMode(model, iComp, deltaLeft, deltaRight, deltaFlankLength, newDIntMax, componentList);
					
		return;
	}
	
	public static double computeMonoInformation(double[] monoBetas, int size) {
		if(size==0) {
			return 0.;
		} else {
			int nChar          = monoBetas.length / size;
			double[] p         = new double[nChar];
			double information = 0;
			double log2        = Math.log(2);
			for(int x=0; x<size; x++) {
				//Computes the probability vector.
				double alphaSum = 0;
				for(int i=0; i<nChar; i++) {
					p[i]      = Math.exp(monoBetas[x*nChar+i]);
					alphaSum += p[i];
				}
				for(int i=0; i<nChar; i++) 
					p[i] /= alphaSum;
				
				//Computes the information.
				information += Math.log(nChar)/log2;
				for(int i=0; i<nChar; i++) {
					if(p[i]>0)
						information += p[i] * Math.log(p[i]) / log2;
				}
			}
			return information;
		}
	}
}
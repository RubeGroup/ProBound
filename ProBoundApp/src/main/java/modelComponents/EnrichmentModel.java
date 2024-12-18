package modelComponents;

import java.util.ArrayList;
import java.util.HashSet;

import org.json.JSONArray;
import org.json.JSONObject;

import proBoundApp.Misc;
import sequenceTools.*;

public abstract class EnrichmentModel extends ModelComponent {
	
	public int nColumns;
	public boolean cumulativeEnrichment;
	String modelType;
	public ArrayList<int[]> iInteractingModes;
	public HashSet<String> modifications;

	//Variables for defining binding modes and interactions.
	ArrayList<BindingMode> allBindingModes;
	ArrayList<BindingModeInteraction> allInteractions;
	public ArrayList<BindingMode> bindingModes;
	public ArrayList<BindingModeInteraction> interactions;
	
	public ArrayList<ArrayList<Integer>> roundSpecificBindingModes, roundSpecificInteractions;
	public boolean noRSBM, noRSint; //indicates if no round-specific binding modes/interactions were specified.
	public boolean[][] roundBindingModeInclusion, roundInteractionInclusion;
	public int nModes = 0, nInteractions = 0;
	public CountTable countTable = null;
	public double concentration;
	

	EnrichmentModel(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModesIn, 
			ArrayList<BindingModeInteraction> allInteractionsIn) {
		super("enrichmentModel");
		
		iComp     = iExpIn;
		
		maxFreezeLevel = 0;
		countTable     = inTable;
		nColumns       = inTable.nColumns;
		
		allBindingModes = new ArrayList<BindingMode>();
		for(BindingMode bm : allBindingModesIn)
			allBindingModes.add(bm);
		
		allInteractions = new ArrayList<BindingModeInteraction>();
		for(BindingModeInteraction bmInt : allInteractionsIn) {
			allInteractions.add(bmInt);
		}

		
		readFromJSON_settings(   config);
		
		
		setFreezingLevel(0);
			
	}
	
	void saveToJSON_settings(JSONObject out) {
		
		String coefficientKey = "modelSettings";
		addEmptyJSON_component_O(out,       coefficientKey, componentKey, iComp);
		JSONObject oEnr   =                  out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oEnr.put("modelType",          modelType);
		
		//Saves the binding mode indices
		JSONArray aBM = new JSONArray();
		for(BindingMode bm : bindingModes)
			aBM.put(bm.iComp);
		oEnr.put("bindingModes", aBM);
		
		//Saves the interaction indices
		JSONArray aInt = new JSONArray();
		for(BindingModeInteraction bmInt : interactions)
			aInt.put(bmInt.iComp);
		oEnr.put("bindingModeInteractions", aInt);
		
		//Saves round-specific binding modes
		JSONArray aRSBM = new JSONArray();
		if(noRSBM)
			aRSBM.put(new JSONArray());
		else {
			for(int iCol=0; iCol<nColumns; iCol++) {
				JSONArray temp = new JSONArray();
				int nRSBMR = roundSpecificBindingModes.get(iCol).size();
				for(int iRSBM=0; iRSBM<nRSBMR; iRSBM++)
					temp.put(bindingModes.get(roundSpecificBindingModes.get(iCol).get(iRSBM)).iComp);
				aRSBM.put(temp);
			}			
		}
		if(aRSBM!=null)
			oEnr.put("roundSpecificBindingModes", aRSBM);
		
		//Saves round-specific binding mode interactions
		JSONArray aRSInt = new JSONArray();
		if(noRSint) //Default empty array
			aRSInt.put(new JSONArray());
		else { //Saves the indices of the binding modes.
			for(int iCol=0; iCol<nColumns; iCol++) {
				JSONArray temp = new JSONArray();
				int nRSIntR = roundSpecificInteractions.get(iCol).size();
				for(int iRSInt=0; iRSInt<nRSIntR; iRSInt++)
					temp.put(interactions.get(roundSpecificInteractions.get(iCol).get(iRSInt)).iComp);
				aRSInt.put(temp);
			}			
		}
		if(aRSInt!=null)
			oEnr.put("roundSpecificBindingModeInteractions", aRSInt);

		
		//Saves the modifications
		JSONArray aMod = new JSONArray();
		for(String mod : modifications)
			aMod.put(mod);
		oEnr.put("modifications", aMod);
		
		//Saves the cumulative enrichment
		oEnr.put("cumulativeEnrichment", cumulativeEnrichment);
		
		//Saves the concentration
		oEnr.put("concentration", concentration);
		
	}
	
	void readFromJSON_settings(JSONObject in) {
		
		String coefficientKey = "modelSettings";
		JSONObject oSettEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		cumulativeEnrichment  = oSettEnr.has("cumulativeEnrichment") ? oSettEnr.getBoolean("cumulativeEnrichment") : true;
		modelType             = oSettEnr.getString("modelType");

		//Building list of binding modes that are included in this enrichment model.
		bindingModes = new ArrayList<BindingMode>();
		JSONArray aBM = oSettEnr.getJSONArray("bindingModes");
		if(aBM.length()==1 && aBM.getInt(0)==-1) //Adds all binding modes.
			for(int iBM=0; iBM<allBindingModes.size(); iBM++)
				bindingModes.add(allBindingModes.get(iBM));
		else { //Adds select binding modes.
			for(int iBM=0; iBM<aBM.length(); iBM++) {
				int iAllBM = aBM.getInt(iBM);
				if(iAllBM < allBindingModes.size() )
					bindingModes.add(allBindingModes.get(iAllBM));
				else
					throw new IllegalArgumentException("Invalid binding mode index index for enrichment model "+iComp+" : " + iAllBM);
			}
		}
		nModes = bindingModes.size();

		//Building list of interactions that are included.
		interactions = new ArrayList<BindingModeInteraction>();
		JSONArray aInt = oSettEnr.getJSONArray("bindingModeInteractions");
		if(aInt.length()==1 && aInt.getInt(0)==-1) //Adds all interactions
			for(int iInt=0; iInt<allInteractions.size(); iInt++)
				interactions.add(allInteractions.get(iInt));
		else { //Adds select binding modes.
			for(int iInt=0; iInt<aInt.length(); iInt++) {
				int iAllInt = aInt.getInt(iInt);
				if(iAllInt < allInteractions.size() )
					interactions.add(allInteractions.get(iAllInt));
				else
					throw new IllegalArgumentException("Invalid binding mode interaction index index for enrichment model "+iComp+" : " + iAllInt);
			}
		}
		nInteractions = interactions.size();
		
		//Builds list of interactions indices:
		iInteractingModes = new ArrayList<int[]>();
		for(int iInt=0; iInt<interactions.size(); iInt++) {
			//Adding new object
			iInteractingModes.add(new int[2]);
			iInteractingModes.get(iInt)[0] = -1;
			iInteractingModes.get(iInt)[1] = -1;
			
			//Identifies the matching interaction.
			for(int iBM=0; iBM<bindingModes.size(); iBM++) {
				
				if(interactions.get(iInt).b0 == bindingModes.get(iBM))
					iInteractingModes.get(iInt)[0] = iBM;
				if(interactions.get(iInt).b1 == bindingModes.get(iBM))
					iInteractingModes.get(iInt)[1] = iBM;
			}
			
			//Throws exception if no matching binding mode was found
			if(iInteractingModes.get(iInt)[0]==-1 || iInteractingModes.get(iInt)[1]==-1) {
				//Creates string with all binding modes in enrichment model
				String bmStr="";
				for(int iBM=0; iBM<bindingModes.size(); iBM++) {
					if(bmStr.length()>0)
						bmStr = bmStr+",";
					bmStr=bmStr+bindingModes.get(iBM).iComp;
					
				}
					
				throw new IllegalArgumentException("Binding mode interaction "+aInt.getInt(iInt)+" in enrichment model "+iComp+" cannot be constructed since the constituent binding modes ("+interactions.get(iInt).b0.iComp+","+interactions.get(iInt).b0.iComp+") are not included in the enrichment model ("+bmStr+").");
			}
		}
		
		//Reads round-specific binding modes
		roundSpecificBindingModes = new ArrayList<ArrayList<Integer>>();
		roundBindingModeInclusion = new boolean[nColumns][nModes];
		JSONArray aRSBM = oSettEnr.has("roundSpecificBindingModes") ? oSettEnr.getJSONArray("roundSpecificBindingModes") : null;
		
		if(aRSBM==null || (aRSBM.length()==1 && aRSBM.getJSONArray(0).length()==0)) {
			noRSBM=true;
			 //Adds all binding modes.
			for(int iR=0; iR<nColumns; iR++) {
				ArrayList<Integer> temp = new ArrayList<Integer>();
				for(int iBM=0; iBM<nModes; iBM++) {
					temp.add(iBM);
					roundBindingModeInclusion[iR][iBM] = true;
				}
				roundSpecificBindingModes.add(temp);
			}
		} else {
			noRSBM=false;
			//Reading the round-specific binding modes.
			//Checks so the number of columns is correct.
			if(aRSBM.length() != nColumns) 
				throw new IllegalArgumentException("The number of columns does not match the length of 'roundSpecificBindingModes' for enrichment model "+iComp+".");
			//Loops over the columns
			for(int iR=0; iR<nColumns; iR++) {
				JSONArray aRSBMR = aRSBM.getJSONArray(iR);
				ArrayList<Integer> temp = new ArrayList<Integer>();
				//Searches for a match for each binding mode
				for(int iBM=0; iBM<aRSBMR.length(); iBM++) {
					int iMatch = -1;
					for(int iTest=0; iTest<nModes; iTest++)
						if (bindingModes.get(iTest).iComp == aRSBMR.getInt(iBM))
							iMatch = iTest;
					if(iMatch==-1) 
						throw new IllegalArgumentException("No matching binding mode found when parsing 'roundSpecificBindingModes' for enrichment model "+iComp+".");
					temp.add(iMatch);
					roundBindingModeInclusion[iR][iMatch] = true;
				}
				roundSpecificBindingModes.add(temp);
			}
		}
		
		//Reads interactions
		roundSpecificInteractions = new ArrayList<ArrayList<Integer>>();
		roundInteractionInclusion = new boolean[nColumns][nInteractions];
		JSONArray aRSInt = oSettEnr.has("roundSpecificBindingModeInteractions") ? oSettEnr.getJSONArray("roundSpecificBindingModeInteractions") : null;
		if(aRSInt==null || (aRSInt.length()==1 && aRSInt.getJSONArray(0).length()==0)) {
			noRSint=true;
			 //Adds all interactions
			for(int iR=0; iR<nColumns; iR++) {
				ArrayList<Integer> temp = new ArrayList<Integer>();
				for(int iInt=0; iInt<nInteractions; iInt++) {
					temp.add(iInt);
					roundInteractionInclusion[iR][iInt] = true;
				}
				roundSpecificInteractions.add(temp);
			}
		} else {
			noRSint=false;
			//Reading the round-specific interactions.
			//Checks so the number of columns is correct.
			if(aRSInt.length() != nColumns) 
				throw new IllegalArgumentException("The number of columns does not match the length of 'roundSpecificBindingModeInteractions' for enrichment model "+iComp+".");
			//Loops over the columns
			for(int iR=0; iR<nColumns; iR++) {
				JSONArray aRSIntR = aRSInt.getJSONArray(iR);
				ArrayList<Integer> temp = new ArrayList<Integer>();
				//Searches for a match for each binding mode interaction
				for(int iInt=0; iInt<aRSIntR.length(); iInt++) {
					int iMatch = -1;
					for(int iTest=0; iTest<nInteractions; iTest++) {
						if (interactions.get(iTest).iComp == aRSIntR.getInt(iInt))
							iMatch = iTest;
					}
					if(iMatch==-1) 
						throw new IllegalArgumentException("No matching binding mode interaction found when parsing 'roundSpecificBindingModeInteractions' for enrichment model "+iComp+".");
					temp.add(iMatch);
					roundInteractionInclusion[iR][iMatch] = true;
				}
				roundSpecificInteractions.add(temp);
			}
		}		
		

				
		//Loads modifications.
		JSONArray aMod = oSettEnr.getJSONArray("modifications");
		modifications  = new HashSet<String>();
		for(int i=0; i<aMod.length(); i++) {
			modifications.add(aMod.getString(i));
		}
		
		//Reads the concentration from the settings JSON object.
		concentration = oSettEnr.has("concentration") ? oSettEnr.getDouble("concentration") : 1;
		
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
	}
	
	@Override
	public void indicateModelChoice(String bestKey) {
		return;
	}
	
	public static EnrichmentModel buildEnrichmentModel(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModes, 
			ArrayList<BindingModeInteraction> allInteractions) {
		String modelType = config.getJSONObject("modelSettings").getJSONArray("enrichmentModel").getJSONObject(iExpIn).getString("modelType");

		if(modelType.equals("SELEX"))
			return new SELEXModel(              config, iExpIn, inTable, allBindingModes, allInteractions);
		else if(modelType.equals("RhoGamma"))
			return new RhoGammaModel(           config, iExpIn, inTable, allBindingModes, allInteractions);
		else if(modelType.equals("Exponential"))
			return new ExponentialModel(        config, iExpIn, inTable, allBindingModes, allInteractions);
		else if(modelType.equals("ExponentialKinetics"))
			return new ExponentialKineticsModel(config, iExpIn, inTable, allBindingModes, allInteractions);
		else if(modelType.equals("InputFreeFraction"))
			return new InputFreeFraction(       config, iExpIn, inTable, allBindingModes, allInteractions);
		else
            throw new IllegalArgumentException("Invalid enrichment model type: " + modelType);
		//return null;		
	}
	
	abstract public void updateDeltaKappaRI(double[] deltaKappaRI, double[] alphaRI, int nRounds);
	abstract public void updateNablaF1(double[] nablaF1, double[] nablaN, double[] alphaRI, int nRounds);
	
	abstract public enrichmentGradient newGradient();
	
	abstract public class enrichmentGradient {
		abstract void addEnrichmentGradientContribution(double[] nablaN, double[] alphaRI);
		abstract void saveToJSON(JSONObject outObject, String coefficientKey);
		abstract void saveToJSON_squared(JSONObject outObject, String coefficientKey);
	};

	public void computeAlphas(ArrayList<SlidingWindow> sw, double[] alphaSeq, double[] alphaInt, double[] alphaRI, 
			ArrayList<ArrayList<ArrayList<Double>>> longAlphaList, LongSequence probe) {

		//Computes the sliding window for all binding modes.
		for(int iBm=0; iBm<nModes; iBm++){	
			ArrayList<ArrayList<Double>> longAlphaListTemp;
			BindingMode oBM = bindingModes.get(iBm);
			if(oBM.includeComponent) {
				if(oBM.k>0){
					if(!bindingModes.get(iBm).swIncludeDi) {
						longAlphaListTemp = sw.get(iBm).slidePN(probe,                       1 );
					} else {
						longAlphaListTemp = sw.get(iBm).slidePN(probe, Math.min(oBM.dInt+1, 2) );
					}

					//Sets alphas on reverse strand to zero if sinlgeStrand is used
					if(oBM.singleStrand) {
						ArrayList<Double> a0 = longAlphaListTemp.get(0), a1 = longAlphaListTemp.get(1);
						double tempSum = 0;
						for (int x = 0; x < a0.size(); x++){
							a1.set(x, 0.);
							tempSum += a0.get(x) ;						
						}
						longAlphaListTemp.get(3).set(0,tempSum);
					}

					//Updates alphas to include the biases.
					double tempSum = 0;
					if(oBM.usePositionBias){
						ArrayList<Double> a0 = longAlphaListTemp.get(0),                 a1 = longAlphaListTemp.get(1);
						double[]          b0 = oBM.positionBiasAlphas.get(iComp).get(0), b1 = oBM.positionBiasAlphas.get(iComp).get(1);

						for (int x = 0; x<b0.length; x++) {
							tempSum += a0.get(x)*b0[x] + a1.get(x)*b1[x];
						}
						
					} else {
						tempSum = longAlphaListTemp.get(3).get(0);
					}

					alphaSeq[iBm] = tempSum;
				} else {
					longAlphaListTemp = null;
					alphaSeq[iBm]     = 1;
				}
			} else {
				longAlphaListTemp = null;
				alphaSeq[iBm]     = 0;
			}


			longAlphaList.set(iBm, longAlphaListTemp);
		}

		//Loops over interactions, computes the interaction-weighted sum
		for(int iInt=0; iInt<nInteractions; iInt++){
			BindingModeInteraction oInt = interactions.get(iInt);
			
			if(oInt.includeComponent) {

				int nF0                                 = oInt.b0.maxFrames.get(iComp)/2;
				int nF1                                 = oInt.b1.maxFrames.get(iComp)/2;
				ArrayList<ArrayList<Double>> alphaList0 = longAlphaList.get(iInteractingModes.get(iInt)[0]);
				ArrayList<ArrayList<Double>> alphaList1 = longAlphaList.get(iInteractingModes.get(iInt)[1]);
				double alphaIntTemp = 0;

				for(int s0=0; s0<2; s0++) {
					for(int s1=0; s1<2; s1++) {
						double[][] intTemp           = oInt.interactionAlphas.get(iComp).get(s0).get(s1);
						ArrayList<Double> alpha0Temp = alphaList0.get(s0);
						ArrayList<Double> alpha1Temp = alphaList1.get(s1);
						for(int x0=0; x0<nF0; x0++) {
							for(int x1=0; x1<nF1; x1++) {
								alphaIntTemp += intTemp[x0][x1] * alpha0Temp.get(x0) * alpha1Temp.get(x1) ;
							}
						}
					}
				}
				alphaInt[iInt] = alphaIntTemp;
				
			} else {
				alphaInt[iInt] = 0;
			}
			
		}

		//Calculates activity weighted alpha_ri
		for(int r=0; r<nColumns; r++){
			alphaRI[r] = 0;

			//Adds contribution from single-binding mode term
			for(Integer iBM : roundSpecificBindingModes.get(r)) {
				BindingMode oBM = bindingModes.get(iBM);
				if(oBM.includeComponent)
					alphaRI[r] += alphaSeq[iBM] * oBM.activityAlphas.get(iComp)[r] * concentration;

			}

//			for(int iBM=0; iBM<nModes; iBM++) {
//				BindingMode oBM = bindingModes.get(iBM);
//				if(oBM.includeComponent)
//					alphaRI[r] += alphaSeq[iBM] * oBM.activityAlphas.get(iComp)[r] * concentration;
//			}

			//Adds contribution from interaction term.
			for(Integer iInt : roundSpecificInteractions.get(r)) {
				BindingModeInteraction oInt = interactions.get(iInt); 
				if(oInt.includeComponent)
					alphaRI[r] += alphaInt[iInt]  * oInt.activityAlphas.get(iComp)[r] * concentration;
			}
				
//			for(int iInt=0; iInt<nInteractions; iInt++) {
//				BindingModeInteraction oInt = interactions.get(iInt); 
//				if(oInt.includeComponent)
//					alphaRI[r] += alphaInt[iInt]  * oInt.activityAlphas.get(iComp)[r] * concentration;
//			}
		}

		if(nModes ==0)
			for(int r=0;r<nColumns; r++)
				alphaRI[r] = 1;
		
	}
	
	public void standardActivityAdjustment(double weight) {
		// Shifts the activities into h so that expected partition function="weight" in each column;

		
				// Special Cases (future):
				// 1. experimentSpecificActivity=TRUE 
				//		DEFAULT: Adjust the activities in the relevant rounds.
				// 2. experimentSpecificActivity=False:
				//		DEFAULT: Adjust the activity with the same shift across all experiments
				//		EXCEPTION:	If some other SELEX model already has includeComponent=True and bindingSaturation=true,
				//		then the activities are already fixed.
						

				double[] expectedPartitionFunction = countTable.computeExpectedPartitionFunction();
				
				System.out.println("Expected partition function before adjustment: "+Misc.formatVector_d(expectedPartitionFunction));

				
				for(int iCol=0; iCol<expectedPartitionFunction.length; iCol++) {

					//Shift required to make the expectedPartitionFunction=weight
					double factor =          weight/expectedPartitionFunction[iCol];
					double shift  = Math.log(factor);
					
					System.out.println("Shift in SELEXModel.activationAdjustment = Log["+weight+"/"+expectedPartitionFunction[iCol]+"] = "+shift);
					
					//Shifts the h.
					countTable.h[iCol]                               -= shift;
					
					//Shifts the binding mode activities.
					for(BindingMode oBM : bindingModes) {
						if(oBM.includeComponent) {
							if(oBM.fitLogActivity)
								oBM.activityBetas.get( iComp)[iCol]  += shift;
							else
								oBM.activityAlphas.get(iComp)[iCol]  *= factor;
						}
					}
					
					//Shifts the binding mode interactions
					for(BindingModeInteraction oInt : interactions) {
						if(oInt.includeComponent) {
							if(oInt.fitLogActivity)
								oInt.activityBetas.get( iComp)[iCol] += shift;
							else
								oInt.activityAlphas.get(iComp)[iCol] *= factor;
						}
					}
				}
				
				for(BindingMode oBM : bindingModes) 
					if(oBM.fitLogActivity)
						oBM.updateAlphas();
				
				
				for(BindingModeInteraction oInt : interactions) 
					if(oInt.fitLogActivity)
						oInt.updateAlphas();


				System.out.println("Expected partition function after adjustment: "+Misc.formatVector_d(countTable.computeExpectedPartitionFunction()));

	}

}

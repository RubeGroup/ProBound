package modelComponents;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.concurrent.*;

import org.json.*;

import base.Array;
import base.MersenneTwisterFast;
import modelComponents.MultiRoundData;
import proBoundApp.Misc;
import sequenceTools.*;

public class CountTable extends ModelComponent  {
	
	//GOAL: Holds the data, computes the likelihood and its gradient.  
	
	//Variables defining dataset.
	private LongSequence.SequenceClass sc;
	public String letterComplement;
	public String letterOrder;
	int nMono, nDi;
	public MultiRoundData fullTable;
	public int l, nColumns;
	public int nReads;
	String countTableFile, inputFileType, rightFlank, leftFlank;
	
	//Variables for testing data
	public MultiRoundData fullTableTest = null;
	public int nReadsTest;
//	String countTableFileTest = null;
	public int[] testFolds = null;
	ArrayList<Integer> stachedProbeIndices = null;
	
	MersenneTwisterFast randGenerator;
	//Variables defining the batch table.
	public ArrayList<Integer> batchProbeIndices;
	public ArrayList<int[]> batchProbeCounts;
	public int[] batchCountPerRound;
	public boolean batchFixedLibrarySize;
	public MultiRoundData batchFullTable;
	ArrayList<int[]> readIndices = null;

	public boolean probeBias;
	public String errorModel;
	boolean usePoisson, useNormal, useLogNormal;
	
	//Round that should be included in the likelihood.
	public int[] modeledColumns;
	
	//Threading information
	//private int iStart, iEnd, nData;
	private int nThreads;
	private int[][] threadRange;
	private ExecutorService pool;

	//Variables defining the enrichment model. 
	public double[] eta, h;

	//Variables for computing the function and its gradients.
	public double weight;
	double lambda;
	public double functionValue; //Function value
	public JSONObject gradient;  //Combined gradient
	public double functionValueSquared; //Sum over the squared value
	public boolean computeVariance; 
	
	private ArrayList<String> trIn, trOut;
	
	//Object defining the enrichment between the columns
	EnrichmentModel enr;
	
	//public boolean verbose = true;

	
	///////////////////////
	// Setting up object //
	///////////////////////
	
	
	//Default constructor (Constructor that doesn't read the data).
	public CountTable(JSONObject config, int iExpIn, LongSequence.SequenceClass scIn, String letterComplementIn, String letterOrderIn, boolean loadData) {
		super("countTable");

		iComp          = iExpIn;
		componentName  = "Count table "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");

		sc               = scIn;
		letterComplement = letterComplementIn;
		letterOrder      = letterOrderIn;
		nMono            = letterOrder.length();
		nDi              = nMono * nMono;
		
		maxFreezeLevel = 0;
		readFromJSON_settings(   config);
		readFromJSON_constraints(config);

		seed_component(config);

		//Sets up threading.
		if(config.has("optimizerSetting") && config.getJSONObject("optimizerSetting").has("nThreads"))
			setNThreads(config.getJSONObject("optimizerSetting").getInt("nThreads"));
		else
			setNThreads(1);

		randGenerator = new MersenneTwisterFast();
		batchFixedLibrarySize = false;
		setFreezingLevel(0);
				
		variationName       = null;
		variationsOptimized = true;
		computeVariance     = false;
		
		if(loadData) {

			//Loads the Training Data
			//TODO: Clean up
			MultiRoundData tempTrainTable = Misc.readDataFile(config, sc, iComp);
//			MultiRoundData tempTrainTable = Misc.readDataFile(config, sc, iComp, false);
			MultiRoundData tempTestTable  = null;
			
//			if(testFolds!=null&&countTableFileTest!=null)
//				throw new IllegalArgumentException("For count table, "+iComp+", both 'countTableFileTest' and 'testFolds' are specified, but at can one should be used.");

			
			//Alternative 1: Splits the training/test data
			if(testFolds!=null&&testFolds.length>0) {
				if(countTableFile==null||!inputFileType.equals("tsv.gz")||countTableFile.length()<6||!countTableFile.endsWith("tsv.gz"))
					throw new IllegalArgumentException("The argument 'testFolds' is only available for 'tsv.gz' input count table files.");
				
				ArrayList<MultiRoundData> splitTables = Misc.splitTableTrainTest(countTableFile.replaceAll("tsv\\.gz$", "folds.gz"), tempTrainTable, testFolds);
				tempTrainTable = splitTables.get(0);
				tempTestTable  = splitTables.get(1);				
			}
			
			//Alternative 2: Loads separate testing data. 
//			if(countTableFileTest!=null)
//				tempTestTable = Misc.readDataFile(config, sc, iComp, true);

			//Loads the training table
			loadDataTable(tempTrainTable, false);
				
			//Loads the testing table
			if(tempTestTable!=null)
				loadDataTable(tempTestTable, true);
			
			nextBatchAll();

		}
			
	}
		
	public void loadDataTable(MultiRoundData dataIn, boolean isTest) {
		
		//Reads data, builds a naive batch.
		if(dataIn.longProbes==null)
			dataIn.encodeLongToLongSequence(sc);

		//Transliterates characters
		if(trIn.size()>0) 
			MultiRoundData.transliterate(dataIn, sc, trIn, trOut);

		//Counts the number of reads
		int nReadsTemp = 0;
		for(int r: modeledColumns)
			nReadsTemp += (int) dataIn.countPerRound[r];
		
		//Checks so the length of the probes are correct
		if(dataIn.longProbes.get(0).getLength() != l)
			throw new IllegalArgumentException("For count table (test="+isTest+")"+iComp+", modelSettings.countTable.variableRegionLength="+l+" does not match the length of the first sequence in the count table ("+dataIn.longProbes.get(0)+").");


		if(isTest) {
			//Saves training data
			fullTableTest      = dataIn;
			nReadsTest         = nReadsTemp;

		} else {
			//Saves training data
			fullTable          = dataIn;
			nReads             = nReadsTemp;
		}
	}
	

	
	//Sets up threading.
	public void setNThreads(int nThreads) {
		this.nThreads = nThreads;
		pool          = Executors.newFixedThreadPool(nThreads);
		
	}
	
	
	public void setEnrichmentModel(EnrichmentModel enrIn) {
		enr           = enrIn;
	}
	
	//////////////////////////////////////////////
	// Function for implementing ModelComponent //
	//////////////////////////////////////////////
		
	public void allocateParameters() {
		
		h     = new double[nColumns];
		eta   = exp_d(h);
		updateAlphas();
		return;
		
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
	
	public void updateAlphas() {
		
		if(eta==null)
			eta =new double[nColumns];
		else
			for(int r=0; r<nColumns; r++)
				eta[r] = 0;
		
		for(int r: modeledColumns) 
			eta[r] = Math.exp(h[r]);
			
		return;
	}

	@Override
	public Hashtable<String,JSONObject> getVariations(JSONObject currentModel, ArrayList<ModelComponent> componentList) {
		
		if(variationsOptimized)
			return null;
		
		//Initial model shift
		/////////////////
		Hashtable<String,JSONObject> newVariations = new Hashtable<String,JSONObject>();
		String testLabel = "initial";
		if(!variationDescription.containsKey(testLabel)) {
			variationDescription.put(testLabel, "Initial model.");
			newVariations.put(testLabel, currentModel);
		}
		if(newVariations.size()>0)
			return newVariations;
		
		//No variation possible.
		variationsOptimized = true;
		return null;

	}
	
	@Override
	public void indicateModelChoice(String bestKey) {
		return;
	}

	@Override
	public void seed_component(JSONObject config) {

		String coefficientKey = "modelSeeding";
		h=null;
		if(config.has(coefficientKey)) {
			JSONObject oSeed = config.getJSONObject(coefficientKey); 
			if( oSeed.has(componentKey) ) {
				JSONObject oEnr = oSeed.getJSONArray(componentKey).getJSONObject(iComp);
				if(oEnr.has("h"))
					h     = readFromJSON_d(oEnr.getJSONArray("h"));
			}
		}
		
		if(h==null)
			h = new double[nColumns];

		updateAlphas();


	}
	
	@Override
	public void activationAdjustment(double weight) {

	}
	
	public double[] computeExpectedPartitionFunction() {
		
		if(!includeComponent)
			return null;
		
		double[] alphaSum = new double[nColumns];
		
		//Sums over binding modes.
		for(int iBM=0; iBM<enr.bindingModes.size(); iBM++) {
			BindingMode oBM = enr.bindingModes.get(iBM);
			
			//Only sums over included binding modes
			if(!oBM.includeComponent)
				continue;
			
			//Expectation value of the scoring matrix.
			double alphaSeq = oBM.computeExpectedAlpha();
			
			//Computes the sum of the position bias/interactions.
			double alphaPBSum = (oBM.usePositionBias && oBM.k>0) ? tr_Ad(oBM.positionBiasAlphas.get(iComp)) : oBM.maxFrames.get(iComp);
			
			//Adds contribution to the alpha sum
			for(int iCol=0; iCol<nColumns; iCol++) {
				for(Integer iRSBM : enr.roundSpecificBindingModes.get(iCol)) //Use correct binding modes in each column
					if(iBM==iRSBM) 
						alphaSum[iCol] += alphaSeq * alphaPBSum * oBM.activityAlphas.get(iComp)[iCol] * enr.concentration;
			}
					
		}
		
		//Sums over interactions
		for(int iInt=0; iInt<enr.interactions.size(); iInt++) {
			BindingModeInteraction oInt = enr.interactions.get(iInt);
			
			//Only sums over included interactions
			if(!oInt.includeComponent)
				continue;
			
			//Expectation value of the scoring matrix.
			double alphaSeq = oInt.b0.computeExpectedAlpha() * oInt.b0.computeExpectedAlpha();

			//Computes the sum of the position bias/interactions.
			double alphaIntSum = tr_AAdd(oInt.interactionAlphas.get(iComp));

			//Adds contribution to the alpha sum
			for(int iCol=0; iCol<nColumns; iCol++)
				for(Integer iRSInt : enr.roundSpecificInteractions.get(iCol))
					if(iInt==iRSInt)
						alphaSum[iCol] += alphaSeq * alphaIntSum * oInt.activityAlphas.get(iComp)[iCol] * enr.concentration;


		}
		
		return alphaSum;
	}
	
	@Override
	int packModel_component(JSONObject packing, int iFirst) {
		
		String coefficientKey = "packing";
		addEmptyJSON_component_O(packing, coefficientKey, componentKey, iComp);

		int iCurr = iFirst;
		
		if(fitComponent) {

			int[] hRange = ModelComponent.constant_i(h, -1);
			for(int i=0; i<modeledColumns.length; i++) 
				hRange[modeledColumns[i]] = iCurr+i+1;
			iCurr += modeledColumns.length;

			saveToJSON_h_i(packing, coefficientKey, hRange);
		}
		return iCurr;
	}

	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_h_d(in,     coefficientKey, zero_d(h));
		
		return;
	}

	@Override
	void saveToJSON_settings(JSONObject out) {

		String coefficientKey = "modelSettings";
		addEmptyJSON_component_O(out,       coefficientKey, componentKey, iComp);
		JSONObject oCT   =                  out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		
		oCT.put("countTableFile",       countTableFile);
		oCT.put("inputFileType",        inputFileType);
		oCT.put("nColumns",             nColumns);
		oCT.put("variableRegionLength", l);
		oCT.put("rightFlank",           rightFlank);
		oCT.put("leftFlank",            leftFlank);
		oCT.put("modeledColumns",       ModelComponent.JSONArrayConvert_i(modeledColumns));
		oCT.put("probeBias",            probeBias);
		oCT.put("errorModel",           errorModel);
	
		//TODO: Clean up
//		if(countTableFileTest!=null)
//			oCT.put("countTableFileTest",       countTableFileTest);
		
		if(testFolds!=null)
			oCT.put("testFolds",        ModelComponent.JSONArrayConvert_i(testFolds));
			
		JSONObject oTR = new JSONObject();
		oCT.put("transliterate", oTR);
		oTR.put("in",            trIn);
		oTR.put("out",           trOut);

	}
	
	@Override
	void readFromJSON_settings(JSONObject in) {

		String coefficientKey = "modelSettings";
		JSONObject oCT        = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		countTableFile        = oCT.getString("countTableFile");
		
		inputFileType         = oCT.getString("inputFileType");
		nColumns              = oCT.getInt("nColumns");
		l                     = oCT.getInt("variableRegionLength");
		rightFlank            = oCT.getString("rightFlank");
		leftFlank             = oCT.getString("leftFlank");
		
		//TODO: Clean up
//		countTableFileTest    = oCT.has("countTableFileTest") ? oCT.getString("countTableFileTest")           : null;
		testFolds             = oCT.has("testFolds")          ? readFromJSON_i(oCT.getJSONArray("testFolds")) : null;

		
		
		probeBias              = oCT.has("probeBias") ? oCT.getBoolean("probeBias") : true;
		errorModel             = oCT.has("errorModel")? oCT.getString("errorModel") : "Poisson";
		
		if(errorModel.equals("Poisson"))
			usePoisson=true;
		if(errorModel.equals("Normal"))
			useNormal=true;
		if(errorModel.equals("LogNormal"))
			useLogNormal=true;
		
		//By default, all columns are modeled.
		modeledColumns = new int[nColumns];
		for(int c=0; c<nColumns; c++)
			modeledColumns[c] = c;
		//Reads the modeled columns from the file if they exist. 
		if(oCT.has("modeledColumns")) {
			JSONArray mc = oCT.getJSONArray("modeledColumns");
			if( !(mc.length()==1&&mc.getInt(0)==-1) )
				modeledColumns = ModelComponent.readFromJSON_i(mc);
		}
		
		trIn  = new ArrayList<String>();
		trOut = new ArrayList<String>();
		if(oCT.has("transliterate")) {
			JSONObject oTR = oCT.getJSONObject("transliterate");
			if(oTR.has("in") && oTR.has("out")) {
				//Gets the in and out JSON arrays
				JSONArray oIn  = (oTR.get("in")  instanceof JSONArray) ? oTR.getJSONArray("in")  : (new JSONArray()).put(oTR.getString("in"));
				JSONArray oOut = (oTR.get("out") instanceof JSONArray) ? oTR.getJSONArray("out") : (new JSONArray()).put(oTR.getString("out"));

				//Checks so 'in' and 'out' have the same number of elements.
				if(oIn.length() != oOut.length())
					throw new IllegalArgumentException("For modelSettings.countTable.transliterate, the length of 'in' does not match the length of 'out.");
				//Loops over elements
				for(int iStr=0; iStr<oIn.length(); iStr++) {
					//Gets the strings
					String sIn  = oIn.getString(iStr);
					String sOut = oOut.getString(iStr);
					//Checks so the in and out strings have the same length 
					if(sIn.length() != sOut.length())
						throw new IllegalArgumentException("For modelSettings.countTable.transliterate, the length of 'in' string "+sIn+" does not match the length of 'out string "+sOut+".");
					//Saves the strings.
					trIn.add(sIn);
					trOut.add(sOut);
				}
			}
		}
	}
	
	@Override
	void saveToJSON_constraints(JSONObject out) {
		
		String coefficientKey = "modelFittingConstraints";
		
		//Creates coefficients if required. 
		if(!out.has(coefficientKey))
			out.append(coefficientKey, new JSONObject());
		JSONObject oCoeff = out.getJSONObject(coefficientKey);
		
		//Creates binding modes if required. 
		if(!oCoeff.has(componentKey))
			oCoeff.append(componentKey, new JSONObject());
		//JSONObject oEnr = oCoeff.getJSONObject(componentKey);
		

	}

	@Override
	void readFromJSON_constraints(JSONObject in) {
		
	}
	
	@Override
	public void regularizationToJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_h_d(in,     coefficientKey, h);

	}
	
	@Override
	void saveToJSON_parameters(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_h_d(in,     coefficientKey, h);

	}

	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {
		JSONObject oEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		h               = readFromJSON_d(oEnr.getJSONArray("h"));
		updateAlphas();
		
	}
	
	//Methods for saving double objects
	void saveToJSON_h_d(JSONObject out, String coefficientKey, double[] h) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("h", JSONArrayConvert_d(h));

		return;
	}
	

	//Methods for saving long objects
	void saveToJSON_h_i(JSONObject out, String coefficientKey, int[] h) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("h", JSONArrayConvert_i(h));

		return;
	}
	
	////////////////////////////////////////////////////////
	// Function for computing function value and gradient //
	////////////////////////////////////////////////////////
	
	//Creates a new batch table given a list of read indices for each column.
	public void nextBatch(ArrayList<int[]> randomReads, double expectedReadCount) {
		
		//TODO: Implement This.
		if(!usePoisson)
			throw new IllegalArgumentException("Batched likelihood evaluator has only been developed for Poisson likelihood and needs to be exteneded to cover the Normal and Log-Normal error functions .");

		
		
		//FIRST TIME ONLY: Creates a list of reads: readIndices[iColumn][iRead] = <index of probe>
		if(readIndices==null) {
			readIndices = new ArrayList<int[]>();
			for(int iCol=0; iCol<nColumns; iCol++) {
				int[] newIndices = new int[fullTable.countPerRound[iCol]];
				readIndices.add(newIndices);
				int iRead = 0;
				for(int iProbe=0; iProbe<fullTable.longProbes.size(); iProbe++) {
					int readCount = fullTable.countTable.get(iProbe)[iCol];
					if(readCount>0) {
						for(int a=0; a<readCount; a++) {
							newIndices[iRead] = iProbe;
							iRead++;
						}
					}
				}
			}
		}
		
		//Variables for the batch table.
		batchProbeIndices  = new ArrayList<Integer>();
		batchProbeCounts   = new ArrayList<int[]>();
		batchCountPerRound = new int[nColumns];
		batchFullTable     = fullTable;
		
		//Map from index of probe in the full table to index of the probe in the batch table.
		HashMap<Integer,Integer> probeToBatchRow  = new HashMap<Integer,Integer>();

		for(int iCol : modeledColumns) {
			//Adds the randomly selected reads to the batch count table.
			
			for(int ir: randomReads.get(iCol)) {
				int iProbe = readIndices.get(iCol)[ir]; //Index of the probe corresponding to the current read.
				
				//If the probe corresponding to the current read does not exist in the batch table, create a new entry.  
				if(!probeToBatchRow.containsKey(iProbe)) {
					batchProbeIndices.add(iProbe);
					batchProbeCounts.add(new int[nColumns]);
					probeToBatchRow.put(iProbe, batchProbeCounts.size()-1);
				}
				
				//Increments the probe count
				batchProbeCounts.get(probeToBatchRow.get(iProbe))[iCol] += 1;
				batchCountPerRound[iCol] ++;
			}
		}
		
		weight = 1. / expectedReadCount;
		
		threadSchedule(nThreads);	
	}
	
	//Creates a new batch given a list of read indices (not grouped by column)
	public void nextBatch(int[] randomUnsortedReads, double expectedReadCount) {
		
		//TODO: Implement This.
		if(!usePoisson)
			throw new IllegalArgumentException("Batched likelihood evaluator has only been developed for Poisson likelihood and needs to be exteneded to cover the Normal and Log-Normal error functions .");

		
		int[] cpr = fullTable.countPerRound;
		
		//A temporary list storing the read indices in each column.  
		ArrayList<ArrayList<Integer>> tempColReads = new ArrayList<ArrayList<Integer>>();
		for(int iCol=0; iCol<nColumns; iCol++)
			tempColReads.add(new ArrayList<Integer>());
		
		//Step 2: Sorts the random reads by columns
		for(int iRead: randomUnsortedReads) {
			int iInCol = iRead;
			//for(int iCol=0; iCol<nColumns; iCol++) {
			for(int iCol : modeledColumns) {
				if(iInCol < cpr[iCol]) {
					tempColReads.get(iCol).add(iInCol);
					break;
				} else {
					iInCol -= cpr[iCol];
				}
			}
		}
		
		//Step 3: Converts the array list to a int[];
		//Creates a random list of the read indices for each column  
		ArrayList<int[]> sortedReads = new ArrayList<int[]>();
		for(int iCol=0; iCol<nColumns; iCol++)  {
			ArrayList<Integer> aI = tempColReads.get(iCol);
			int[] tableCol = new int[aI.size()];
			for(int iI=0; iI<aI.size(); iI++)
				tableCol[iI] = aI.get(iI);
			sortedReads.add(tableCol);
		}

		nextBatch(sortedReads, expectedReadCount);
	}
	
	//Creates a new batch using all reads.
	public void nextBatchAll() {
		//Uses all reads in the table
		batchFullTable     = fullTable;
		batchProbeCounts   = batchFullTable.countTable;
		batchCountPerRound = batchFullTable.countPerRound;
		if(stachedProbeIndices==null) {
			batchProbeIndices  = new ArrayList<Integer>();
			for(int i=0; i<batchProbeCounts.size(); i++)
				batchProbeIndices.add(i);
		} else {
			batchProbeIndices   = stachedProbeIndices;
			stachedProbeIndices = null;
		}
		
		//Schedules the threads.
		threadSchedule(nThreads);
		
		if(usePoisson)
			weight = 1.0 / nReads;
		else 
			weight = 1.0 / batchFullTable.longProbes.size();

	}
	
	//Creates a new batch using all testing reads.
	public void nextBatchAllTest() {
		//Uses all reads in the table
		batchFullTable      = fullTableTest;
		batchProbeCounts    = batchFullTable.countTable;
		batchCountPerRound  = batchFullTable.countPerRound;
		stachedProbeIndices = batchProbeIndices;
		batchProbeIndices   = new ArrayList<Integer>();
		for(int i=0; i<batchProbeCounts.size(); i++)
			batchProbeIndices.add(i);
		
		//Schedules the threads.
		threadSchedule(nThreads);
		
		if(usePoisson)
			weight = 1.0 / nReadsTest;
		else 
			weight = 1.0 / batchFullTable.longProbes.size();

	}
	
	//Creates a new batch given the total number of desired reads.
	public void nextBatch(int nDataIn) {
		
		//TODO: Implement This.
		if(!usePoisson)
			throw new IllegalArgumentException("Batched likelihood evaluator has only been developed for Poisson likelihood and needs to be exteneded to cover the Normal and Log-Normal error functions .");

		if(fullTable.countTable!=null) {

			//Counts the number of reads
			nReads             = 0;
			for(int r: modeledColumns)
				nReads += (int) fullTable.countPerRound[r];
			
			//Determines the whole dataset should be used.
			int nData;
			if(nDataIn>=nReads) //Makes sure that the batch is not lager than the number of probe sequences
				nData = -1;
			else if(nDataIn==0) //Makes sure we use at least one probe.
				nData = 1;
			else
				nData = nDataIn;
			
			if(nData == -1) {
				nextBatchAll();
				return;

			} else {
				
				//int nReadsTot = (int) Array.sum(fullTable.countPerRound);
				if(batchFixedLibrarySize) {
					//CASE 1: Number of reads per column proportional to sequencing depth 
					//////////////////////////////////////////////////////////////////////
					//Counts the total number of reads
					int[] tempCountPerRound  = new int[nColumns];
					//Creates a random list of the read indices for each column  
					ArrayList<int[]> randomReads = new ArrayList<int[]>();
					int nReadsCum = 0; //Cumulative number of reads in the columns
					for(int iCol=0; iCol<nColumns; iCol++) {
						//Only pull reads from the modeled columns.
						if(Array.memberQ(modeledColumns, iCol)) {

							//Step 1: Determines how many reads we should have in each column. 
							int lastNReadsCum        = nReadsCum;
							nReadsCum               += fullTable.countPerRound[iCol];
							
							tempCountPerRound[iCol]  =  ((int) Math.floor(1.0*nData*nReadsCum    /nReads)) 
										              - ((int) Math.floor(1.0*nData*lastNReadsCum/nReads));
							//Step 2: Creates a list of read indices.
							randomReads.add(Misc.randomSample(fullTable.countPerRound[iCol], tempCountPerRound[iCol], randGenerator));
						} else {
							//Add no reads if the columns isn't modeled.
							randomReads.add(new int[0]);
						}
					}
					
					//Creates count table.
					nextBatch(randomReads, nData);
					return;

				} else {
					//CASE 2: Reads drawn randomly regardless of round 
					//////////////////////////////////////////////////
					nextBatch(Misc.randomSample(nReads, nData, randGenerator), nData);
					return;
				}
			}
			
			
		} else {
			nReads        = 0;
			weight        = 0;
			
			//Schedules the threads.
			threadSchedule(nThreads);	
			return;
		}

		
	}
	
	public void writeBatchTable(String outTable) {
		
		//Loops over probes.
		PrintStream originalStream = System.out;
		PrintStream fileStream     = null;
		try {
			fileStream = new PrintStream(new FileOutputStream(outTable, false));
			System.setOut(fileStream);
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create batch table file: "+outTable);
			e.printStackTrace();
			System.exit(1);
		}
		
		for(int iProbe=0; iProbe<batchProbeCounts.size(); iProbe++)
			System.out.println(batchFullTable.longProbes.get(batchProbeIndices.get(iProbe)).toString()+"\t"+Misc.formatVector_i(batchProbeCounts.get(iProbe), "\t", "", ""));
			
		//Returns to original stream.
		System.setOut(originalStream);
		fileStream.close();
				
	}

	// Sets up threading
	////////////////////
	public void threadPoolShutdown() {
		pool.shutdown();
		while (!pool.isShutdown()) {}
		return;
	}
	
	private void threadSchedule(int nThreads) {
		
		//Threading setup for cycling the dats
		threadRange		                = new int[nThreads][2];
		int nProbes                     = batchProbeCounts.size();
		
		for(int iThread=0; iThread<nThreads; iThread++ ) {
			threadRange[iThread][0]     = ((int) Math.floor((1.0*(iThread  )/nThreads)*nProbes )); 
			threadRange[iThread][1]     = ((int) Math.floor((1.0*(iThread+1)/nThreads)*nProbes )); 
			
			// Checks if the thread range contains zero reads: threadRange[iThread][0]=threadRange[iThread][1]
			// could imply that the the thread contains zero probes or that it contains all probes. 
			if( nProbes<nThreads && threadRange[iThread][0]==threadRange[iThread][1] ) {
				threadRange[iThread][0] = -1;
				threadRange[iThread][1] = -1;
			}
		}
		
	}		

	// Functions for launch threads for computing function value and the gradient
	/////////////////////////////////////////////////////////////////////////////
	
	//Updates the value.
	public void updateValue() throws InterruptedException, ExecutionException {
		
		functionValue               = 0;
		if(computeVariance)
			functionValueSquared    = 0;  
		List<Future<JSONObject>> threadOutput;
		Set<Callable<JSONObject>> tasks = new HashSet<Callable<JSONObject>>(nThreads);
		JSONObject threadResult;
		
		//Assign Threads
		for (int i=0; i<nThreads; i++) 
			if(threadRange[i][0]!= -1)
				tasks.add(new ThreadedFunctionEvaluator(threadRange[i][0], threadRange[i][1]));

		//Launch Threads
		threadOutput = pool.invokeAll(tasks);
		//Sum up value and return
		for (Future<JSONObject> currentThread : threadOutput) {
			threadResult   = currentThread.get();
			//Note: The sign flips: threadResults.functionValue is log-likelihood, 'functionValue' is negative log-likelihood  
			functionValue -= threadResult.getDouble("functionValue");
			if(computeVariance) {
				functionValueSquared += threadResult.getDouble("functionValueSquared");
			}
		}
	}

	//Updates the gradient.
	public void updateGradient() throws InterruptedException, ExecutionException {
		
		functionValue = 0.0;
		gradient = createEmptyJSONModel("gradient");
		
		//Sums up the squared values.
		if(computeVariance) {
			functionValueSquared    = 0;
			ModelComponent.addEmptyJSONObject(gradient, "gradientSquared");
		}

		//NOTE: the gradient object contains "functionValue"
		
		List<Future<JSONObject>> threadOutput;
		Set<Callable<JSONObject>> tasks = new HashSet<Callable<JSONObject>>(nThreads);
		JSONObject threadResult;
		
		//Assign Threads
		for (int i=0; i<nThreads; i++) 
			if(threadRange[i][0]!= -1)
				tasks.add(new ThreadedGradientEvaluator(threadRange[i][0], threadRange[i][1]));

			//Launch Threads
		threadOutput = pool.invokeAll(tasks);
		//Sum up value and return
		for (Future<JSONObject> currentThread : threadOutput) {
			threadResult              = currentThread.get();
			//Note: The sign flips: threadResults.functionValue is log-likelihood, 'functionValue' is negative log-likelihood (same for gradient).
			functionValue            -= threadResult.getDouble("functionValue");
			ModelComponent.addJSONObject_O(gradient.getJSONObject("gradient"), threadResult.getJSONObject("gradient"), -1);
			
			//Sums up the squared values.
			if(computeVariance) {
				ModelComponent.addJSONObject_O(gradient.getJSONObject("gradientSquared"), threadResult.getJSONObject("gradientSquared"), 1);
				functionValueSquared += threadResult.getDouble("functionValueSquared");
			}
		}
	}
	
	private void computePRI(double[] pRI, double[] kappaRI, double[] deltaKappaRI, double[] etaKappa) {
		
		//Calculates Kappa
		if(enr.cumulativeEnrichment){
			double runningKappaProduct = 1.0;
			//TODO: Should the product start at r=1? In theory this the first factor should cancel out in PRI, but numerical errors might decrease. 
			for(int r=0; r<nColumns; r++){
				if(r>0)
					runningKappaProduct *= deltaKappaRI[r];
				kappaRI[r] = runningKappaProduct;
			}
		} else {
			for(int r=0; r<nColumns; r++){
				kappaRI[r] = deltaKappaRI[r];
			}
		}

		//Updates pRI using only the columns specified in modeledColumns.	
		double zProbe = 0.0;
		for(int r=0; r<nColumns; r++)
			pRI[r] = 0;
		
		for(int r: modeledColumns) {
			etaKappa[r] = eta[r] * kappaRI[r];
			zProbe += etaKappa[r];
		}
		
		for(int r: modeledColumns) 
			pRI[r] = etaKappa[r] / zProbe;

	}
	
	public void writePCTable(String predictedTable) /*throws Exception*/ {
		
		int nRounds = fullTable.countPerRound.length;

		double[] alphaSeq             = new double[enr.nModes];
		double[] alphaInt             = new double[enr.nInteractions];
		double[] alphaRI              = new double[nRounds];
		ArrayList<ArrayList<ArrayList<Double>>> longAlphaList = new ArrayList<ArrayList<ArrayList<Double>>>(); 
		for(int bm=0; bm<enr.nModes; bm++)
			longAlphaList.add(null);

		double[] deltaKappaRI         = new double[nRounds];
		double[] kappaRI              = new double[nRounds];
		double[] pRI                  = new double[nRounds];
		double[] etaKappa             = new double[nRounds];
		
		ArrayList<SlidingWindow> sw  = new ArrayList<SlidingWindow>();
		for(int bm=0; bm<enr.nModes; bm++)
			sw.add(enr.bindingModes.get(bm).getSlidingWindow(iComp, enr.modifications));

		//Loops over probes.
		PrintStream original	= System.out;
		try {
			System.setOut(new PrintStream(new FileOutputStream(predictedTable, false)));
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create predicted-table file at this location: "+predictedTable);
			e.printStackTrace();
			System.exit(1);
		}

		int nDataPoints		= fullTable.longProbes.size();
		
		//Loops over all probes (ignoring possible batching and multithreading.)
		for(int iProbe=0; iProbe<nDataPoints; iProbe++) {
			//Updates alpha values.
			enr.computeAlphas(sw, alphaSeq, alphaInt, alphaRI, longAlphaList, fullTable.longProbes.get(iProbe));
			//Calculates deltaKappa:
			enr.updateDeltaKappaRI(deltaKappaRI, alphaRI, nRounds);
			//Computes pRI.
			computePRI(pRI, kappaRI, deltaKappaRI, etaKappa);
			
			//Prefroms reverse transliteration to get the original sequence.
			String probeSeq = fullTable.longProbes.get(iProbe).toString();
			for(int iTr=0; iTr<trIn.size(); iTr++) 
				probeSeq = probeSeq.replaceAll(trOut.get(iTr), trIn.get(iTr));
			
			//Multiplies pRI by the number of reads.
			/*System.out.println(probeSeq + "\t" + Misc.formatVectorE_d(Array.scalarMultiply(pRI, nRow),"\t","","", 10));*/
			
			
			
			////
			if(usePoisson || useNormal) {
								
				if(probeBias) {
					
					//Counts the reads in the modeled columns.
					long nRow = 0;
					for(int r:modeledColumns)
						nRow += fullTable.countTable.get(iProbe)[r];

					System.out.println(probeSeq + "\t" + Misc.formatVectorE_d(Array.scalarMultiply(pRI, nRow),"\t","","", 10));
				} else {
					
					System.out.println(probeSeq + "\t" + Misc.formatVectorE_d(etaKappa,"\t","","", 10));

				}
			} else {
								
				if(probeBias) {
					
					double geometricMeanN = 1.0, geometricmeanEtaKappa = 1.0;
					for(int r: modeledColumns) {
						geometricMeanN        *= fullTable.countTable.get(iProbe)[r];
						geometricmeanEtaKappa *= etaKappa[r];
					}
					geometricMeanN        = Math.pow(geometricMeanN,        1/modeledColumns.length);
					geometricmeanEtaKappa = Math.pow(geometricmeanEtaKappa, 1/modeledColumns.length);


					System.out.println(probeSeq + "\t" + Misc.formatVectorE_d(Array.scalarMultiply(pRI, geometricMeanN/geometricmeanEtaKappa),"\t","","", 10));

				} else {

					System.out.println(probeSeq + "\t" + Misc.formatVectorE_d(etaKappa,"\t","","", 10));
					
				}
			}
			
			///
		} 
		
		//Returns to original stream.
		System.setOut(original);
		
		return;
	}
	
	
	public void writeAlphaTable(String outAlphaFile, boolean writeBMInt) /*throws Exception*/ {
		
		int nRounds                   = fullTable.countPerRound.length;
		double[] alphaSeq             = new double[enr.nModes];
		double[] alphaInt             = new double[enr.nInteractions];
		double[] alphaRI              = new double[nRounds];
		
		ArrayList<ArrayList<ArrayList<Double>>> longAlphaList
		                              = new ArrayList<ArrayList<ArrayList<Double>>>();
		
		for(int bm=0; bm<enr.nModes; bm++)
			longAlphaList.add(null);

		ArrayList<SlidingWindow> sw   = new ArrayList<SlidingWindow>();
		for(int bm=0; bm<enr.nModes; bm++)
			sw.add(enr.bindingModes.get(bm).getSlidingWindow(iComp, enr.modifications));

		//Loops over probes.
		int nSequences		          = fullTable.longProbes.size();
		PrintStream original          = System.out;
		
		try {
			System.setOut(new PrintStream(new FileOutputStream(outAlphaFile, false)));
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create output table at this location: "+outAlphaFile);
			e.printStackTrace();
			System.exit(1);
		}
		
		//Loops over sequences, computes affinity sums (alphas), and writes to file
		for(int iProbe=0; iProbe<nSequences; iProbe++) {
			LongSequence currentProbeSeq = fullTable.longProbes.get(iProbe);
			enr.computeAlphas(sw, alphaSeq, alphaInt, alphaRI, longAlphaList, currentProbeSeq);

			//Writes the table
			String outString   = fullTable.longProbes.get(iProbe).toString();
			if(writeBMInt) {
				//Writes the alpha sum for each binding mode and interaction separately
				
				//Adds contribution from single-binding mode term
				double[] ouAlphas = new double[nColumns];
				for(int iBM=0; iBM<enr.nModes; iBM++) {
					BindingMode oBM = enr.bindingModes.get(iBM);
					for(int r=0; r<nColumns; r++){
						if(oBM.includeComponent) {
							ouAlphas[r] = alphaSeq[iBM] * oBM.activityAlphas.get(iComp)[r] * enr.concentration;
						} else {
							ouAlphas[r] = 0;
						}
					}
					outString += "\t" + Misc.formatVector_d(ouAlphas, ",", "", "", 8);
				}

				//Adds contribution from interaction term.
				for(int iInt=0; iInt<enr.nInteractions; iInt++) {
					BindingModeInteraction oInt = enr.interactions.get(iInt); 
					for(int r=0; r<nColumns; r++){
						if(oInt.includeComponent)
							ouAlphas[r] =  alphaInt[iInt] * oInt.activityAlphas.get(iComp)[r] * enr.concentration;
						else
							alphaRI[r]  = 0;
					}
					outString += "\t" + Misc.formatVector_d(ouAlphas, ",", "", "", 8);
				}

				
			} else {
				//Only writes the full alpha sum
				outString += "\t" + Misc.formatVector_d(alphaRI, "\t", "", "", 8);
			}

			System.out.println(outString);
		}
		
		//Returns to original stream.
		System.setOut(original);
		
		return;
	}
	
	public void writeBindingModeAlphas(String outAlphaFile, String outFormat) {
		
		//Creates a list of relevant sliding windows and binding modes
		ArrayList<BindingMode> bmList     = new ArrayList<BindingMode>();
		ArrayList<SlidingWindow> swList   = new ArrayList<SlidingWindow>();
		for(BindingMode oBM: enr.bindingModes) {
			if(oBM.includeComponent && oBM.k>0) {
				bmList.add(oBM);
				swList.add(oBM.getSlidingWindow(iComp, enr.modifications));
			}
		}
		
		//Loops over probes.
		int nSequences		          = fullTable.longProbes.size();
		PrintStream original          = System.out;
		
		//Opens output file
		if(!outAlphaFile.equals("-")) {
			try {
				System.setOut(new PrintStream(new FileOutputStream(outAlphaFile, false)));
			} catch (FileNotFoundException e) {
				System.out.println("Cannot create output table at this location: "+outAlphaFile);
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		
		//Loops over sequences, computes affinity sums (alphas), and writes to file
		ArrayList<ArrayList<Double>> longAlphaList;
		for(int iProbe=0; iProbe<nSequences; iProbe++) {
			
			LongSequence currentProbeSeq = fullTable.longProbes.get(iProbe);
			
			System.out.print(currentProbeSeq.toString());
			
			for(int iBM=0; iBM<bmList.size(); iBM++) {
				
				//Scores the sequence
				BindingMode oBM   = bmList.get(iBM);
				SlidingWindow oSW = swList.get(iBM);
				boolean ss        = oBM.singleStrand;
				
				if(!oBM.swIncludeDi) {
					longAlphaList = oSW.slidePN(currentProbeSeq,                       1 );
				} else {
					longAlphaList = oSW.slidePN(currentProbeSeq, Math.min(oBM.dInt+1, 2) );
				}
				String colString = "";
				int nDigits      = 5;
				if(outFormat.equals("max")) {
					double maxValue  = Array.max(longAlphaList.get(0));
					if(!ss)
						maxValue = Math.max(maxValue,  Array.max(longAlphaList.get(1)));
					colString        = String.format("%."+nDigits+"e", maxValue);

				} else if(outFormat.equals("sum")) {
					double sumValue  = Array.sum(longAlphaList.get(0));
					if(!ss)
							sumValue+= Array.sum(longAlphaList.get(1));
					colString        = String.format("%."+nDigits+"e", sumValue);
					
				} else if(outFormat.equals("mean")) {
					int n = longAlphaList.get(0).size();
					double meanValue =                 (Array.sum(longAlphaList.get(0))/n);
					if(!ss) 
						meanValue = 0.5 * (meanValue + (Array.sum(longAlphaList.get(1))/n));
					colString        = String.format("%."+nDigits+"e", meanValue);

				} else if(outFormat.equals("profile")) {
					Collections.reverse(longAlphaList.get(1));
					colString        = "" +   Misc.formatVectorE_d(                      longAlphaList.get(0),         ",", "", "", 5);
					if(!ss)
						colString   += "\t" + Misc.formatVectorE_d(                      longAlphaList.get(1),         ",", "", "", 5);
					else
						colString   += "\t" + Misc.formatVectorE_d(ModelComponent.zero_d(longAlphaList.get(1).size()), ",", "", "", 5);
					
					
					
				} else {
					System.err.println("ERROR: Invalid output format: '"+outFormat+"'.");
					System.exit(1);
				}
				System.out.print("\t"+colString);
			}
			System.out.println("");
		}
		
		//Returns to original stream.
		System.setOut(original);
		
		return;
		
	}
	
	//Uses the count table to compute a weight:
	// (maxStrandMean / minStrandMean) * maxValue
	public HashMap<Integer,Double> getEmpericalBindingModeActivities() {
		
		//TODO: Consider updating to respect round-specifc binding modes/interactions.
		HashMap<Integer,Double> out = new HashMap<Integer,Double>();
				
		//Loops over binding modes (that are relevant for this count table)
		ArrayList<ArrayList<Double>> longAlphaList;
		//BindingMode oBM = enr.bindingModes.get(iBMInt);
		for(BindingMode oBM : enr.bindingModes) {

			//Checks so the binding mode is included and has size>0, skips if not.
			if(!oBM.includeComponent || oBM.k<=1)
				continue;

			//Gets the sliding window.
			int tempFL        = oBM.flankLength;
			oBM.flankLength   = 0;
			SlidingWindow oSW = oBM.getSlidingWindow(iComp, enr.modifications);
			oBM.flankLength   = tempFL;

			//Loops over probes
			int countSum = 0;
			double forwardSum=0, reverseSum=0, maxVal=0, activitySum = 0;
			
			double[] vSumFwrd = new double[nColumns], vSumRvrs = new double[nColumns];

			
			for(int iProbe=0; iProbe<fullTable.longProbes.size(); iProbe++) {
				LongSequence currentProbeSeq = fullTable.longProbes.get(iProbe);
				
				//Scores using sliding window (that doesn't go into the flanks).
				if(!oBM.swIncludeDi) {
					longAlphaList = oSW.slidePN(currentProbeSeq,                       1 );
				} else {
					longAlphaList = oSW.slidePN(currentProbeSeq, Math.min(oBM.dInt+1, 2) );
				}

				//Computes the max value and the strand sums
				double probeMeanForward = Array.mean(longAlphaList.get(0));
				double probeMeanReverse = Array.mean(longAlphaList.get(1));
				double weight           = 0;
				boolean bmSelected      = false;
				for(int iCol : modeledColumns) {
					//Gets the counts, computes a count-based probe weight 
					int cnt         = fullTable.countTable.get(iProbe)[iCol];
					weight         += cnt;
					
					//Counts the reads and sums the count-weighted activities, ignoring the first 
					//round if cumulativeEnrichment=true:  
					if( !(enr.cumulativeEnrichment&&iCol==0) ) {
						activitySum += cnt * oBM.activityAlphas.get(iComp)[iCol];
						countSum    += cnt;
						bmSelected   = true;
					}
					
					
				}
				
				//For computing the round-wise ratio
				for(int iCol=0; iCol<nColumns; iCol++) {
					vSumFwrd[iCol] += probeMeanForward * fullTable.countTable.get(iProbe)[iCol];
					vSumRvrs[iCol] += probeMeanReverse * fullTable.countTable.get(iProbe)[iCol];
				}
				
				//Computes the appropriate round-and-activity weighted sums
				forwardSum += weight * probeMeanForward;
				reverseSum += weight * probeMeanReverse;
				
				//Computes maximum value (only keeping probes with non-zero counts in a modeled column)
				if(bmSelected) 
						maxVal = Math.max(maxVal, Math.max(Array.max(longAlphaList.get(0)), Array.max(longAlphaList.get(1))));
				
			}
			
			//Ratios for verbose output


			
			//Computes the strand ratio using trend across rounds
			double e1=0, eRatio=0, eRound=0, eRatioRound=0, eRound2=0;
			for(int iCol : modeledColumns) {
				if(fullTable.countPerRound[iCol]>0) {
					double logRatio = Math.log(vSumRvrs[iCol]/vSumFwrd[iCol]);
					
					e1          += fullTable.countPerRound[iCol];
					eRatio      += fullTable.countPerRound[iCol]*logRatio;
					eRound      += fullTable.countPerRound[iCol]*iCol;
					eRatioRound += fullTable.countPerRound[iCol]*iCol*logRatio;
					eRound2     += fullTable.countPerRound[iCol]*iCol*iCol;
				}
			}
			double covRatioRound = (eRatioRound/e1) - (eRatio/e1)*(eRound/e1);
			double varRound      = (eRound2/e1)     - (eRound/e1)*(eRound/e1);
			double slope         = covRatioRound / varRound;
			double modelRatio    = Math.exp(-Math.abs(slope)*(Array.max(modeledColumns)-Array.min(modeledColumns)));

			//Computes the ratio strand ratio by taking read-weighted strand average
			double frwdRevRatio = Math.min(reverseSum/forwardSum, forwardSum/reverseSum);

			//Identifies the most extreme ratio
			double strandRatio = Math.min(modelRatio, frwdRevRatio);
			
			//Deals with None:
			if(Double.isNaN(strandRatio))
				strandRatio = 0;
			
			//Computes the mean position bias (ignoring positions overlapping flanks..
			double meanPB = 1;
			if(oBM.positionBiasAlphas!=null && iComp<oBM.positionBiasAlphas.size()) {
				int fl            = oBM.flankLength;
				double[] pbFwrd   = oBM.positionBiasAlphas.get(iComp).get(0);
				double[] bpRvrs   = oBM.positionBiasAlphas.get(iComp).get(1);
				double pbFwrdMean = Array.mean(Array.copyOfRange(pbFwrd, fl, pbFwrd.length-fl));
				double pbRvrsMean = Array.mean(Array.copyOfRange(bpRvrs, fl, pbFwrd.length-fl));
				meanPB = (pbFwrdMean + pbRvrsMean) / 2;
			}
	
			out.put(oBM.iComp, maxVal * meanPB * (activitySum / countSum) * strandRatio);

			
			if(verbose) {

				double[] vRatio = new double[nColumns];
				for(int iCol=0; iCol<nColumns; iCol++) {
					if(fullTable.countPerRound[iCol]>0) {
						vRatio[iCol] = vSumRvrs[iCol]/vSumFwrd[iCol];
					}
				}
				
				System.out.println("iBM = "+oBM.iComp+", ratio1 = "+Misc.formatVector_d(vRatio)+" => "+modelRatio+", ratio2 = "+frwdRevRatio+" => min="+strandRatio);
				System.out.println("         max="+maxVal+", PB="+(meanPB) + ", act="+(activitySum / countSum)+", ratio="+strandRatio);
			}
		}
				
		return out;
	}
	
	
	public class ThreadedFunctionEvaluator implements Callable<JSONObject>{
		private int startIdx;
		private int endIdx;
		
		public ThreadedFunctionEvaluator(int startIdx, int endIdx) {
			this.startIdx				= startIdx;
			this.endIdx					= endIdx;
		}
		
		
		//@Override
		public JSONObject call() throws Exception {
			
			int nRounds                   = batchCountPerRound.length;

			double[] alphaSeq             = new double[enr.nModes];
			double[] alphaInt             = new double[enr.nInteractions];
			double[] alphaRI              = new double[nRounds];
			ArrayList<ArrayList<ArrayList<Double>>> longAlphaList 
			                              = new ArrayList<ArrayList<ArrayList<Double>>>(); 
			for(int bm=0; bm<enr.nModes; bm++)
				longAlphaList.add(null);

			double[] deltaKappaRI         = new double[nRounds];
			double[] kappaRI              = new double[nRounds];
			double[] pRI                  = new double[nRounds];
			double[] etaKappa             = new double[nRounds];
			
			ArrayList<SlidingWindow> sw = new ArrayList<SlidingWindow>();
			for(int bm=0; bm<enr.nModes; bm++)
				sw.add(enr.bindingModes.get(bm).getSlidingWindow(iComp, enr.modifications));

			double functionValue          = 0.0;
			double functionValueSquared   = 0.0;
			
			//Loops over probes.
			if(startIdx!=-1) {
				for(int i=startIdx; i<endIdx; i++) {
					//Updates alpha values.
					enr.computeAlphas(sw, alphaSeq, alphaInt, alphaRI, longAlphaList, batchFullTable.longProbes.get(batchProbeIndices.get(i)));
												  
					//Calculates deltaKappa:
					enr.updateDeltaKappaRI(deltaKappaRI, alphaRI, nRounds);
					
					//Computes pRI.
					computePRI(pRI, kappaRI, deltaKappaRI, etaKappa);
				
					double probeValue = 0;
					if(usePoisson) {
						if(probeBias) {
							for(int r: modeledColumns)
								probeValue += batchProbeCounts.get(i)[r] * Math.log(pRI[r]);
						} else {
							for(int r: modeledColumns)
								probeValue += batchProbeCounts.get(i)[r] *  Math.log(etaKappa[r]) - etaKappa[r];
						}
					} else if(useNormal) {
						if(probeBias) {
							double nI = 0;
							for(int r: modeledColumns)
								nI += batchProbeCounts.get(i)[r];
							for(int r: modeledColumns)
								probeValue += - Math.pow(batchProbeCounts.get(i)[r] - nI*pRI[r],2);
						} else {
							for(int r: modeledColumns)
								probeValue += - Math.pow(batchProbeCounts.get(i)[r] - etaKappa[r],2);
						}
					} else {
						if(probeBias) {
							double geometricMeanN = 1.0, geometricmeanEtaKappa = 1.0;
							for(int r: modeledColumns) {
								geometricMeanN        *= batchProbeCounts.get(i)[r];
								geometricmeanEtaKappa *= etaKappa[r];
							}
							geometricMeanN        = Math.pow(geometricMeanN,        1/modeledColumns.length);
							geometricmeanEtaKappa = Math.pow(geometricmeanEtaKappa, 1/modeledColumns.length);
							for(int r: modeledColumns)
								probeValue += - Math.pow(Math.log((batchProbeCounts.get(i)[r]/geometricMeanN) / (etaKappa[r]/geometricmeanEtaKappa)),2);
						} else {
							for(int r: modeledColumns)
								probeValue += - Math.pow(Math.log(batchProbeCounts.get(i)[r]/etaKappa[r]),2);
							
						}
					}
					
					functionValue        += probeValue;
					functionValueSquared += probeValue*probeValue;
					



				}
			}
			
			
			JSONObject outObject = new JSONObject();
			outObject.put("functionValue",        functionValue);
			outObject.put("functionValueSquared", functionValueSquared);
			
			return outObject;
		}
	}
	
	public class ThreadedGradientEvaluator implements Callable<JSONObject>{
		private int startIdx;
		private int endIdx;
		private double[] deltaN                               = null;
		private double[] nablaN                               = null;

		//Gradient variables.
		ArrayList<double[]> monoGradients                     = null;
		private ArrayList<ArrayList<double[]>>   diGradients  = null;
		private ArrayList<double[]>   activityGradients       = null;
		private ArrayList<ArrayList<double[]>> positionBiasGradient   
															  = null;
		private ArrayList<ArrayList<ArrayList<double[][]>>> interactionGradients
														      = null;
		private ArrayList<double[]>   intActivityGradients    = null;
		private ArrayList<ArrayList<double[]>> 
		                             monoModGradients         = null;
		private ArrayList<ArrayList<ArrayList<double[]>>>
		                             diModGradients           = null;
		private double[] hGradients	                          = null;
		private EnrichmentModel.enrichmentGradient enrGradient= null;
		
		//Squared-gradient variables.
		ArrayList<double[]> monoGradientsSquared                     = null;
		private ArrayList<ArrayList<double[]>>   diGradientsSquared  = null;
		private ArrayList<double[]>   activityGradientsSquared       = null;
		private ArrayList<ArrayList<double[]>> positionBiasGradientSquared   
															         = null;
		private ArrayList<ArrayList<ArrayList<double[][]>>> interactionGradientsSquared
														             = null;
		private ArrayList<double[]>   intActivityGradientsSquared    = null;
		private ArrayList<ArrayList<double[]>> 
		                             monoModGradientsSquared         = null;
		private ArrayList<ArrayList<ArrayList<double[]>>>
		                             diModGradientsSquared           = null;
		private double[] hGradientsSquared	                         = null;

		public ThreadedGradientEvaluator(int startIdx, int endIdx) {
			this.startIdx				= startIdx;
			this.endIdx					= endIdx;
			
		}
				
		//@Override
		public JSONObject call() throws Exception {
			
			double functionValue	    = 0;
			double functionValueSquared = 0;
			
			ArrayList<BindingMode> bindingModes            = enr.bindingModes;
			ArrayList<BindingModeInteraction> interactions = enr.interactions;
			int nModes                                     = enr.nModes;
			int nInteractions                              = enr.nInteractions;
			
			//Allocates variables for binding mode gradients
			activityGradients          = new ArrayList<double[]>();
			positionBiasGradient       = new ArrayList<ArrayList<double[]>>();
			enrGradient                = enr.newGradient();
			ArrayList<SlidingWindow> sw = new ArrayList<SlidingWindow>();

			ArrayList<double[]> psamGradients        = new ArrayList<double[]>();
			ArrayList<double[]> psamGradientsProbe   = new ArrayList<double[]>();
			ArrayList<double[]> psamGradientsSquared = new ArrayList<double[]>();
			
			
			int nMonoInt = (int) Math.pow(2, sc.getbBits());
			int nDiInt   = nMonoInt * nMonoInt;
			for(int bm=0; bm<nModes; bm++){
				BindingMode oBM = bindingModes.get(bm);
				
				sw.add(oBM.getSlidingWindow(iComp, enr.modifications));
				
				if(oBM.fitMono) {
					//The number of dinucleotide parameters: Sum_{d=1..dInt} (k-d) = dInt*(2*k-dInt-1)/2
					int nParameters = nMonoInt*oBM.k + nDiInt * Math.min(oBM.dIntMax, Math.max(oBM.k-1, 0)) * (2*oBM.k - 1 - Math.min(oBM.dIntMax, Math.max(oBM.k-1, 0)) ) / 2;
					psamGradients.add(       new double[nParameters]);
					psamGradientsSquared.add(new double[nParameters]);
					
					//Only computes the full gradient for each probe if the dinucleotide interactions actually are included
					int nPsamGradCompute       = ( oBM.dInt==0 || (oBM.dInt>0 && oBM.swIncludeDi) ) ? nParameters : nMonoInt*oBM.k;
					psamGradientsProbe.add(  new double[nPsamGradCompute]);
					
				} else {
					psamGradients.add(       null);
					psamGradientsProbe.add(  null);
					psamGradientsSquared.add(null);
				}
				
				if(oBM.fitActivity) {
					activityGradients.add(zero_d(oBM.activityAlphas.get(iComp)));
				} else {
					activityGradients.add(null);
				}
				
				if(oBM.fitPositionBias)
					positionBiasGradient.add(zero_Ad(oBM.positionBiasAlphas.get(iComp)));
				else
					positionBiasGradient.add(null);					
			}
			
			//Allocates variables for interaction gradients
			interactionGradients       = new ArrayList<ArrayList<ArrayList<double[][]>>>();
			intActivityGradients       = new ArrayList<double[]>();
			for(int i=0; i<nInteractions; i++){
				BindingModeInteraction oInt = interactions.get(i);
				if(oInt.fitComponent) {
					interactionGradients.add(zero_AAdd(oInt.interactionAlphas.get(iComp)));
					intActivityGradients.add(zero_d(oInt.activityAlphas.get(iComp)));
				} else {
					interactionGradients.add(null);
					intActivityGradients.add(null);
				}
			}

			//Allocates variables for rho, gamma, and h gradients
			if(fitComponent)
				hGradients               = zero_d(h);
			else
				hGradients               = null;
			
			//Allocates the squared-gradient variables.
			if(computeVariance) {

				if(activityGradients!=null)
					activityGradientsSquared    = ModelComponent.clone_Ad(activityGradients);
				if(positionBiasGradient!=null)
					positionBiasGradientSquared = ModelComponent.clone_AAd(positionBiasGradient);
				if(interactionGradients!=null)
					interactionGradientsSquared = ModelComponent.clone_AAAdd(interactionGradients);
				if(intActivityGradients!=null)
					intActivityGradientsSquared = ModelComponent.clone_Ad(intActivityGradients);
				if(hGradients!=null)
					hGradientsSquared           = ModelComponent.clone_d(hGradients);
			}
			
			
			//Temporary variables
			double[] alphaSeq             = new double[nModes];
			double[] alphaInt             = new double[nInteractions];
			double[] alphaRI              = new double[nColumns];
			double[] deltaKappaRI         = new double[nColumns];
			double[] kappaRI              = new double[nColumns];
			double[] etaKappa             = new double[nColumns];
			double[] pRI                  = new double[nColumns];	
			double[] gradientWeight       = new double[nModes];
			double[] nablaF1              = new double[nColumns];
			

			//Indexing: [bindingMode][0,1,2][position]:
			ArrayList<ArrayList<ArrayList<Double>>> longAlphaList = new ArrayList<ArrayList<ArrayList<Double>>>();			
			for(int bm=0; bm<nModes; bm++)
				longAlphaList.add(null);
			
			//The PSAM gradient has the form W * X, where W=psamGradientAlphas are coefficients and X is an indicator.
			ArrayList<ArrayList<ArrayList<Double>>> psamGradientWeights = new ArrayList<ArrayList<ArrayList<Double>>>();
			for(int bm=0; bm<nModes; bm++) {
				BindingMode oBM = bindingModes.get(bm);
				int maxFrames = oBM.maxFrames.get(iComp);
				
				if(oBM.fitMono && oBM.k > 0) {
					ArrayList<ArrayList<Double>> tempWeights = new ArrayList<ArrayList<Double>>();
					tempWeights.add(new ArrayList<Double>());
					tempWeights.add(new ArrayList<Double>());

					for(int x=0; x<maxFrames; x++){
						tempWeights.get(0).add(0.);
						tempWeights.get(1).add(0.);
					}
					tempWeights.add(new ArrayList<Double>()); //Constant value
					tempWeights.get(2).add(0.);
					psamGradientWeights.add(tempWeights);
				} else {
					psamGradientWeights.add(null);
				}
			}
			
			deltaN         = new double[nColumns];
			nablaN         = new double[nColumns];
			
			  ///////////////////////////////
			 // START OF LOOP OVER PROBES //
			///////////////////////////////
			if(startIdx!=-1) {
				for(int i=startIdx; i<endIdx; i++) {
				//do {
					
					LongSequence currentSeq = batchFullTable.longProbes.get(batchProbeIndices.get(i));

					//Updates alpha values.
					enr.computeAlphas(sw, alphaSeq, alphaInt, alphaRI, longAlphaList, currentSeq);

					//Calculates deltaKappa:
					enr.updateDeltaKappaRI(deltaKappaRI, alphaRI, nColumns);

					//Computes pRI.
					computePRI(pRI, kappaRI, deltaKappaRI, etaKappa);

					//The delta/nabla N weight..
					int[] countRI = batchProbeCounts.get(i);
					Integer nI = 0;
					for(int r: modeledColumns) nI       += countRI[r];                 // nI = number of reads across rounds
					for(int r: modeledColumns) deltaN[r] = 0;                          // Only uses the modeled columns when computing...
					
					double probeValue    = 0.0;

					if(usePoisson) { // Computes log L and d(log L) / dkappa
						if(probeBias) { 
							//Poisson error, With Probe Bias
							for(int r : modeledColumns) {
								
								probeValue += countRI[r] * Math.log(pRI[r]);
								deltaN[r]   = (countRI[r] - nI * pRI[r]);
							}
						} else {
							//Poisson error, No Probe Bias
							for(int r : modeledColumns) {
								probeValue += countRI[r] *  Math.log(etaKappa[r]) - etaKappa[r];
								deltaN[r]   = (countRI[r] - etaKappa[r]);
							}
						}
					} else if(useNormal) {
						if(probeBias) {
							//Normal error, With Probe Bias
							
							//Computes delta n
							double[] dn = new double[nColumns];
							for(int r : modeledColumns)
								dn[r] = countRI[r] - nI*pRI[r];
							//Computes delta-delta n
							double[] ddn = new double[nColumns];
							double sumDnPRI=0;
							for(int r : modeledColumns)
								sumDnPRI += dn[r]*pRI[r];
							for(int r : modeledColumns)
								ddn[r] = dn[r] - sumDnPRI;
							
							for(int r : modeledColumns) {
								probeValue += - Math.pow(countRI[r] - nI*pRI[r],2);
								deltaN[r]  += 2*nI*pRI[r]*ddn[r];
 							}
						} else {
							//Normal error, No Probe Bias
							for(int r: modeledColumns) { 
								probeValue += - Math.pow(countRI[r] - etaKappa[r],2);
								deltaN[r] = 2*(countRI[r] - etaKappa[r])*etaKappa[r];
							}
						}
					} else if (useLogNormal) {
						if(probeBias) { 
							//Log-Normal error, With Probe Bias
							double geometricMeanN = 1.0, geometricmeanEtaKappa = 1.0;
							for(int r: modeledColumns) {
								geometricMeanN        *= countRI[r];
								geometricmeanEtaKappa *= etaKappa[r];
							}
							geometricMeanN        = Math.pow(geometricMeanN,        1/modeledColumns.length);
							geometricmeanEtaKappa = Math.pow(geometricmeanEtaKappa, 1/modeledColumns.length);
							for(int r : modeledColumns) {
								double m = countRI[r]/geometricMeanN;
								double q = etaKappa[r]/geometricmeanEtaKappa;
								probeValue += - Math.pow(Math.log(m/q),2);
								deltaN[r]  += 2*Math.log(m/q);
							}
						} else {
							//Log-Normal error, No Probe Bias
							for(int r : modeledColumns) {
								probeValue += - Math.pow(Math.log(countRI[r]/etaKappa[r]),2);
								deltaN[r]   = 2*Math.log(countRI[r]/etaKappa[r]);
							}
						}
					}
					

					if(enr.cumulativeEnrichment){
						double cumulativeN = 0.0 ;
						for(int r=nColumns-1;r>=0;r--) {
							cumulativeN += deltaN[r];
							nablaN[r] = cumulativeN; // Upper cumulative of deltaN
						}
					} else {
						for(int r=0; r<nColumns; r++) nablaN[r] = deltaN[r];
					}
					
	
					functionValue            += probeValue;
					if(computeVariance)
						functionValueSquared += probeValue*probeValue;
					//Computes nablaN * f^1.
					enr.updateNablaF1(nablaF1, nablaN, alphaRI, nColumns);

					  ////////////////////////////
					 // Binding mode gradients // 
					////////////////////////////
					for(int bm=0; bm<nModes; bm++){

						BindingMode oBM = bindingModes.get(bm);
						if(oBM.fitMono) {

							if(oBM.k > 0) {
								for(int s=0; s<2; s++) {

									ArrayList<Double> alphaTemp = longAlphaList.get(bm).get(s);
									int nF = oBM.maxFrames.get(iComp)/2;

									for(int x=0; x<nF; x++) {
										double tempWeight = 0;

										//First term - single-mode with position bias
										if(oBM.usePositionBias) {
											double[] pbTemp             = oBM.positionBiasAlphas.get(iComp).get(s);
											for(int r=0; r<nColumns; r++) {
												if(enr.roundBindingModeInclusion[r][bm])
													tempWeight += nablaF1[r] * oBM.activityAlphas.get(iComp)[r] * enr.concentration * alphaTemp.get(x) * pbTemp[x];
											}
										} else {
											for(int r=0; r<nColumns; r++) {
												if(enr.roundBindingModeInclusion[r][bm])
													tempWeight += nablaF1[r] * oBM.activityAlphas.get(iComp)[r] * enr.concentration * alphaTemp.get(x);
											}
										}

										//Second term - gradient contribution due to interactions.	
										for(int iI=0; iI<nInteractions; iI++) {
											BindingModeInteraction oInt = interactions.get(iI);

											if(oInt.includeComponent) {

												int bm0 = enr.iInteractingModes.get(iI)[0];
												int bm1 = enr.iInteractingModes.get(iI)[1];
												int nf0 = oInt.b0.maxFrames.get(iComp)/2;
												int nf1 = oInt.b1.maxFrames.get(iComp)/2;

												for(int r=0; r<nColumns; r++){
													if(enr.roundInteractionInclusion[r][iI]) {
														double intActConcTemp = oInt.activityAlphas.get(iComp)[r] * enr.concentration;

														if(bm == bm0) {
															for(int sP=0; sP<2; sP++) {
																ArrayList<Double> alphaTempP = longAlphaList.get(bm1).get(sP);
																double[][] alphaIntTemp      = oInt.interactionAlphas.get(iComp).get(s).get(sP);
																for(int xP=0; xP<nf1; xP++)
																	tempWeight += nablaF1[r] * alphaIntTemp[x][xP] * intActConcTemp * alphaTemp.get(x)  * alphaTempP.get(xP);
															}
														}

														if(bm == bm1) {
															for(int sP=0; sP<2; sP++) {
																ArrayList<Double> alphaTempP = longAlphaList.get(bm0).get(sP);
																double[][] alphaIntTemp      = oInt.interactionAlphas.get(iComp).get(sP).get(s);
																for(int xP=0; xP<nf0; xP++) {
																	tempWeight += nablaF1[r] * alphaIntTemp[xP][x] * intActConcTemp * alphaTempP.get(xP) * alphaTemp.get(x);
																}
															}
														}
													}
												}
											}
										}

										//Sets the gradient weight.
										psamGradientWeights.get(bm).get(s).set(x, tempWeight);

									}
								}
								
								//Actually updates the gradient.
								double[] tempGrad       = psamGradientsProbe.get(bm);
								double[] currGrad       = psamGradients.get(bm);
								double[] currGradSqured = psamGradientsSquared.get(bm);

								//Resets probe gradient
								for(int iTemp=0; iTemp<tempGrad.length; iTemp++)
									tempGrad[iTemp]            = 0;

								//Computes the gradients
								sw.get(bm).swGradient(currentSeq.getValue(), currentSeq.getLength(), 1., psamGradientWeights.get(bm), tempGrad);
								
								//Adds the gradient to the running sum
								for(int iTemp=0; iTemp<tempGrad.length; iTemp++)
									currGrad[iTemp]           += tempGrad[iTemp];
								
								//Saves the squared gradient
								if(computeVariance)
									for(int iTemp=0; iTemp<tempGrad.length; iTemp++)
										currGradSqured[iTemp] += tempGrad[iTemp] * tempGrad[iTemp];
								
							}
						}


						//Calculates the gradient of the activities. Note: This is the activity gradient, not the log(Activity) gradient!!!!
						if(oBM.fitActivity) {
							for(int r=0; r<nColumns; r++)
								if(enr.roundBindingModeInclusion[r][bm])
									activityGradients.get(bm)[r]            +=          nablaF1[r] * alphaSeq[bm] * enr.concentration;
							
							//Computes the variance of the activityygradient.
							if(computeVariance)
								for(int r=0; r<nColumns; r++)
									if(enr.roundBindingModeInclusion[r][bm])
										activityGradientsSquared.get(bm)[r] += Math.pow(nablaF1[r] * alphaSeq[bm] * enr.concentration, 2);
							
						}

						//Calculates gradient of position bias:
						if(oBM.fitPositionBias && oBM.k > 0){

							//  gradientWeight = f^1 * nablaN[r] * a[bm,r]  = dL / dLog[Alpha] 					
							gradientWeight[bm] = 0;
							for(int r=0; r<nColumns; r++)
								if(enr.roundBindingModeInclusion[r][bm])
									gradientWeight[bm] += nablaF1[r] * oBM.activityAlphas.get(iComp)[r] * enr.concentration;

							ArrayList<Double> a0 = longAlphaList.get(bm).get(0),             a1 = longAlphaList.get(bm).get(1);
							double[]          b0 = oBM.positionBiasAlphas.get(iComp).get(0), b1 = oBM.positionBiasAlphas.get(iComp).get(1);
							int al=a0.size();

							double[] pbg0        = positionBiasGradient.get(bm).get(0);
							double[] pbg1        = positionBiasGradient.get(bm).get(1);
							for (int x = 0; x < al; x++) {
								pbg0[x]         += gradientWeight[bm] * a0.get(x) * b0[x];
								pbg1[x]         += gradientWeight[bm] * a1.get(x) * b1[x];
							}
							
							//Computes the variance of the position bias gradient.
							if(computeVariance) {
								
								double[] pbg0_sq = positionBiasGradientSquared.get(bm).get(0);
								double[] pbg1_sq = positionBiasGradientSquared.get(bm).get(1);
								for (int x = 0; x < al; x++) {
									pbg0_sq[x]  += Math.pow(gradientWeight[bm] * a0.get(x) * b0[x], 2);
									pbg1_sq[x]  += Math.pow(gradientWeight[bm] * a1.get(x) * b1[x], 2);
								}
								
							}
						}
					}

					  ////////////////////////////////////////
					 // Binding-mode interaction gradients // 
					////////////////////////////////////////
					for(int iI=0; iI<nInteractions; iI++) {
						BindingModeInteraction oInt = interactions.get(iI);
						if(oInt.fitComponent) {
							int bm1 = enr.iInteractingModes.get(iI)[0];
							int bm2 = enr.iInteractingModes.get(iI)[1];
							int nf1 = oInt.b0.maxFrames.get(iComp)/2;
							int nf2 = oInt.b1.maxFrames.get(iComp)/2;

							//Motif interaction gradient.
							for(int s1=0; s1<2; s1++) {
								for(int s2=0; s2<2; s2++) {
									ArrayList<Double> alpha1                 = longAlphaList.get(bm1).get(s1);
									ArrayList<Double> alpha2                 = longAlphaList.get(bm2).get(s2);
									double[][] alphaIntTemp                  = oInt.interactionAlphas.get(iComp).get(s1).get(s2);
									double[]   alphaIntActTemp               = oInt.activityAlphas.get(iComp);

									//Gradient of interactions 
									double[][] alphaIntGrad                  = interactionGradients.get(iI).get(s1).get(s2);
									double[][] alphaIntGradSquared 	         =computeVariance ? interactionGradientsSquared.get(iI).get(s1).get(s2) : null;
									for(int x1=0; x1<nf1; x1++) {
										for(int x2=0; x2<nf2; x2++) { 
											double gradSum = 0.0;
											for(int r=0; r<nColumns; r++) 
												if(enr.roundInteractionInclusion[r][iI])
													gradSum                     += nablaF1[r] * alphaIntTemp[x1][x2] * alphaIntActTemp[r] * enr.concentration * alpha1.get(x1) * alpha2.get(x2);
											alphaIntGrad[x1][x2]            += gradSum;
											if(computeVariance)
												alphaIntGradSquared[x1][x2] += gradSum * gradSum;
										}
									}

									//Interaction-activity Gradients...
									double[] alphaIntActGrad                 =                   intActivityGradients.get(iI);
									double[] alphaIntActGradSquared          = computeVariance ? intActivityGradientsSquared.get(iI) : null;
									for(int r=0; r<nColumns; r++) {
										if(enr.roundInteractionInclusion[r][iI]) {
											double gradSum = 0.0;
											for(int x1=0; x1<nf1; x1++)
												for(int x2=0; x2<nf2; x2++)
													gradSum                     += nablaF1[r] * alphaIntTemp[x1][x2]                      * enr.concentration* alpha1.get(x1) * alpha2.get(x2);
											alphaIntActGrad[r] += gradSum;
											if(computeVariance)
												alphaIntActGradSquared[r]       += gradSum * gradSum;
										}
									}
								}
							}
						}
					}

					  ///////////////////////////
					 // Count table gradients // 
					///////////////////////////
					if(fitComponent) {

						//Updates h-gradient
						for(int r : modeledColumns)
						//for(int r=0; r<nColumns; r++)
							hGradients[r]            += deltaN[r];
						
						//Computes the squared h-gradient
						if(computeVariance)
							//for(int r=0; r<nColumns; r++)
							for(int r : modeledColumns)
								hGradientsSquared[r] += deltaN[r] * deltaN[r];
					}

					  ////////////////////////////////
					 // Enrichment model gradients // 
					////////////////////////////////
					if(enr.fitComponent) {
						//Adds contribution to the gradient of the enrichment model(such as the rho- and gamma-gradients)
						enrGradient.addEnrichmentGradientContribution(nablaN, alphaRI); //NOTE: This also computes the squared gradient
					}
				}
			}

			  /////////////////////////////
			 // END OF LOOP OVER PROBES //
			/////////////////////////////
			
			//Converts the activity gradient into a log-activity gradient.
			for(int bm=0; bm<bindingModes.size(); bm++){
				BindingMode oBM = bindingModes.get(bm);

				if( oBM.fitActivity && oBM.fitLogActivity ) {
					double[] actAlpha = oBM.activityAlphas.get(iComp);
					for(int iA=0; iA<actAlpha.length; iA++)
						activityGradients.get(bm)[iA]  *= actAlpha[iA];
					
					if(computeVariance)
						for(int iA=0; iA<actAlpha.length; iA++)
							activityGradientsSquared.get(bm)[iA] *= actAlpha[iA] * actAlpha[iA];
				}
			}
			
			for(int iInt=0; iInt<interactions.size(); iInt++) {
				BindingModeInteraction oInt = interactions.get(iInt);
				if( oInt.fitComponent && oInt.fitLogActivity ) {
					double[] intActAlpha = oInt.activityAlphas.get(iComp);

					for(int iA=0; iA<intActAlpha.length; iA++)
						intActivityGradients.get(iInt)[iA]  *= intActAlpha[iA];
					
					if(computeVariance)
						for(int iA=0; iA<intActAlpha.length; iA++)
							intActivityGradientsSquared.get(iInt)[iA] *= intActAlpha[iA] * intActAlpha[iA];
				}
			}


			//Converts flattened PSAM gradient to mono- and di-nucleotide arrays.
			monoGradients              = new ArrayList<double[]>();
			diGradients                = new ArrayList<ArrayList<double[]>>();
			//Creates modification-gradient
			monoModGradients           = new ArrayList<ArrayList<double[]>>();
			diModGradients             = new ArrayList<ArrayList<ArrayList<double[]>>>();


			//Squared mono and di gradients.
			if(computeVariance){
				monoGradientsSquared    = new ArrayList<double[]>();
				diGradientsSquared      = new ArrayList<ArrayList<double[]>>();
				monoModGradientsSquared = new ArrayList<ArrayList<double[]>>();
				diModGradientsSquared   = new ArrayList<ArrayList<ArrayList<double[]>>>();
			}
			
			for(int bm=0; bm<nModes; bm++){
				BindingMode oBM                                     = bindingModes.get(bm);
				double[] tempMono                                   = null;
				ArrayList<double[]> tempDi                          = null;
				ArrayList<double[]> tempModMono_Ad                  = null;
				ArrayList<ArrayList<double[]>> tempModDi_AAd        = null;
				//Variables for squared
				double[] tempMonoSquared                            = null;
				ArrayList<double[]> tempDiSquared                   = null;
				ArrayList<double[]> tempModMonoSquared_Ad           = null;
				ArrayList<ArrayList<double[]>> tempModDiSquared_AAd = null;
				
				if(oBM.fitMono) {
					
					if(sw.get(bm)!=null && psamGradients.get(bm)!=null) {
						psamGradients.set(           bm,sw.get(bm).convertGradient(psamGradients.get(bm),        letterOrder));
						if(computeVariance)
							psamGradientsSquared.set(bm,sw.get(bm).convertGradient(psamGradientsSquared.get(bm), letterOrder));
					}
					
					//Reorganizes mononucleotide gradient contributions
					int iFirst            = 0;
					int iLast             = oBM.monoBetas.length;
					tempMono              = Arrays.copyOfRange(psamGradients.get(bm), iFirst, iLast);
					if(computeVariance)
						tempMonoSquared   = Arrays.copyOfRange(psamGradientsSquared.get(bm), iFirst, iLast);
					iFirst                = iLast;
					
					//Saves mononucleotide modification gradient componentss.
					tempModMono_Ad        = ModelComponent.zero_Ad(oBM.monoBetaShiftValues);
					tempModMonoSquared_Ad = ModelComponent.zero_Ad(oBM.monoBetaShiftValues);
					for(int iMod=0; iMod<oBM.modifications.size(); iMod++) {
						if(enr.modifications.contains(oBM.modifications.get(iMod))) {
							int[] monoIndices                 = oBM.monoBetaShiftIndices.get(iMod);

							//Extracts mononucleotide modifications
							double[] tempModMono_d            = tempModMono_Ad.get(iMod);
							for(int x=0; x<tempModMono_d.length; x++)
								tempModMono_d[x]              = tempMono[monoIndices[x]];
							
							//Extracts squared mononucleotide concentrations
							if(computeVariance) {
								double[] tempModMonoSquared_d = tempModMonoSquared_Ad.get(iMod);
								for(int x=0; x<tempModMonoSquared_d.length; x++)
									tempModMonoSquared_d[x]   = tempMonoSquared[monoIndices[x]];
							}
						}
					}
					

					if(oBM.fitDi) {
						//Reorganizes dinucleotide gradient contributions
						tempDi               = new ArrayList<double[]>();
						tempDiSquared        = new ArrayList<double[]>();
						for(int iD=0; iD<oBM.diBetas.size(); iD++) {
							iLast            = iFirst + oBM.diBetas.get(iD).length;
							tempDi.add(Arrays.copyOfRange(psamGradients.get(bm), iFirst, iLast));
							if(computeVariance)
								tempDiSquared.add(Arrays.copyOfRange(psamGradientsSquared.get(bm), iFirst, iLast));
							iFirst           = iLast;
						}
						
						//Saves dinucleotide modification gradient components.
						tempModDi_AAd        = ModelComponent.zero_AAd(oBM.diBetaShiftValues);
						tempModDiSquared_AAd = ModelComponent.zero_AAd(oBM.diBetaShiftValues);
						for(int iMod=0; iMod<oBM.modifications.size(); iMod++) {
							if(enr.modifications.contains(oBM.modifications.get(iMod))) {
								ArrayList<int[]> diIndices_Ad                 = oBM.diBetaShiftIndices.get(iMod);

								//Extracts variance of dinucleotide modifications
								ArrayList<double[]> tempModDi_Ad              = tempModDi_AAd.get(iMod);
								for(int d=0; d<tempModDi_Ad.size(); d++)
									for(int x=0; x<tempModDi_Ad.get(d).length; x++) 
										tempModDi_Ad.get(d)[x]                = tempDi.get(d)[diIndices_Ad.get(d)[x]];
								
								//Extracts squared variance of dinucleotide modifications 
								if(computeVariance) {
									ArrayList<double[]> tempModDiSquared_Ad   = tempModDiSquared_AAd.get(iMod);
									for(int d=0; d<tempModDiSquared_Ad.size(); d++)
										for(int x=0; x<tempModDiSquared_Ad.get(d).length; x++) 
											tempModDiSquared_Ad.get(d)[x]     = tempDiSquared.get(d)[diIndices_Ad.get(d)[x]];

								}
							}
						}
					}
				}
				
				monoGradients.add(          tempMono);
				diGradients.add(            tempDi);
				monoModGradients.add(       tempModMono_Ad);
				diModGradients.add(         tempModDi_AAd);
				
				if(computeVariance) {
					monoGradientsSquared.add(   tempMonoSquared);
					diGradientsSquared.add(     tempDiSquared);
					monoModGradientsSquared.add(tempModMonoSquared_Ad);
					diModGradientsSquared.add(  tempModDiSquared_AAd);
				}
				

			}
			
			// Combines the gradients into a single JSON object.
			String gradientKey = "gradient";
			String squaredKey  = "gradientSquared";

			JSONObject outObject = createEmptyJSONModel(gradientKey);
			if(computeVariance)
				ModelComponent.addEmptyJSONObject(outObject, squaredKey);

			//Binding modes
			for(int iBM=0; iBM<nModes; iBM++) {
				BindingMode oBM = bindingModes.get(iBM);
				if( oBM.fitComponent ) {
					//Savse the binding mode gradient
					oBM.addZeroJSON_component(                                   outObject, gradientKey);
					if( oBM.fitMono )         oBM.saveToJSON_mono_d(             outObject, gradientKey, monoGradients.get(iBM));
					if( oBM.fitDi )           oBM.saveToJSON_di_Ad(              outObject, gradientKey, diGradients.get(iBM));
					if( oBM.fitPositionBias ) oBM.saveToJSON_positionBias_Ad(    outObject, gradientKey, positionBiasGradient.get(iBM), iComp);
					if( oBM.fitActivity )     oBM.saveToJSON_activity_d(         outObject, gradientKey, activityGradients.get(iBM), iComp);

					oBM.saveToJSON_modification_Ad_AAd (outObject, gradientKey, 
						oBM.fitMono ? monoModGradients.get(iBM) : null, 
						oBM.fitDi   ? diModGradients.get(iBM)   : null);
					
					//Saves the squared gradient  
					if(computeVariance) {
						oBM.addZeroJSON_component(                               outObject, squaredKey);
						if( oBM.fitMono )         oBM.saveToJSON_mono_d(         outObject, squaredKey, monoGradientsSquared.get(iBM));
						if( oBM.fitDi )           oBM.saveToJSON_di_Ad(          outObject, squaredKey, diGradientsSquared.get(iBM));
						if( oBM.fitPositionBias ) oBM.saveToJSON_positionBias_Ad(outObject, squaredKey, positionBiasGradientSquared.get(iBM), iComp);
						if( oBM.fitActivity )     oBM.saveToJSON_activity_d(     outObject, squaredKey, activityGradientsSquared.get(iBM), iComp);

						oBM.saveToJSON_modification_Ad_AAd (                     outObject, squaredKey, 
							oBM.fitMono ? monoModGradientsSquared.get(iBM) : null, 
							oBM.fitDi   ? diModGradientsSquared.get(iBM)   : null);
					}
				}
			}
			
			//Binding mode interactions.
			for(int iInt=0; iInt<nInteractions; iInt++) {
				BindingModeInteraction oInt = interactions.get(iInt);
				if( oInt.fitComponent) {
					//Saves the interaction gradient
					oInt.addZeroJSON_component(          outObject, gradientKey);
					oInt.saveToJSON_interaction_AAdd(    outObject, gradientKey, interactionGradients.get(iInt), iComp);
					oInt.saveToJSON_activity_d(          outObject, gradientKey, intActivityGradients.get(iInt), iComp);
					
					//Saves the squared gradient
					if(computeVariance) {
						oInt.addZeroJSON_component(      outObject, squaredKey);
						oInt.saveToJSON_interaction_AAdd(outObject, squaredKey,  interactionGradientsSquared.get(iInt), iComp);
						oInt.saveToJSON_activity_d(      outObject, squaredKey,  intActivityGradientsSquared.get(iInt), iComp);
					}
				}
			}
			
			//Count table parameters
			if(fitComponent) {
				//Saves the h-gradient
				addZeroJSON_component(                   outObject, gradientKey);
				saveToJSON_h_d(                          outObject, gradientKey, hGradients);
			}
			
			//Enrichment model parameters
			if(enr.fitComponent) {
				//Gradient from the current enrichment model. 
				enr.addZeroJSON_component(               outObject, gradientKey);
				enrGradient.saveToJSON(                  outObject, gradientKey);
				
				//Saves the squared gradient
				if(computeVariance) {
					enr.addZeroJSON_component(           outObject, squaredKey);
					enrGradient.saveToJSON_squared(      outObject, squaredKey);
				}
			}
			outObject.put("functionValue",               functionValue);
			//Saves the squared function value
			if(computeVariance)
				outObject.put("functionValueSquared",    functionValueSquared);
			
			return outObject;

		}
	}
	
	static public void printJSONObjectCoefficients(JSONObject model, String coefficientKey, int iExp) {
		
		JSONObject oTable = model.getJSONObject(coefficientKey).getJSONArray("countTable").getJSONObject(iExp);
		System.out.println("h:                 "+Misc.formatVector_d(ModelComponent.readFromJSON_d(oTable.getJSONArray("h"))));
		
	}
}
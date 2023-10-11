package proBoundApp;

import configBuilder.BuildingFunction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.*;

import modelOptimization.*;
import sequenceTools.LongSequence;
import sequenceTools.SlidingWindow;

public  class App 
{
	
	
	
	static String proBoundDir       = null;
    public static String generalSchemaFile = null;
    public static String configSchemaFile  = null;
    public static String builderSchemaFile = null;
    
	static JSONModel config;
	
	public static void main(String[] args) {
				
		String version = "v1.0.0";
		//System.out.println("> Starting ProBound " +version);
        
		if(System.getenv("PROBOUND_DIR")==null) {
			System.err.println("ERROR: The environmental variable PROBOUND_DIR has not been set.");
			System.exit(1);
		}
		String proBoundDir = System.getenv("PROBOUND_DIR");
		generalSchemaFile  = proBoundDir + "/config/" + "schema.general.json";
		configSchemaFile   = proBoundDir + "/config/" + "schema.config.json";
		builderSchemaFile  = proBoundDir + "/config/" + "schema.builder.json";
		
		//  *********************
		// ** PARSING ARGUMENTS ** 
		//  *********************
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();

		Option configOption = new Option("c", "config", true, "config file path");
		configOption.setRequired(true);
		options.addOption(configOption);
		
		Option builderOption = new Option("b", "builder", false, "the config");
		builderOption.setRequired(false);
		options.addOption(builderOption);
		
		HelpFormatter formatter  = new HelpFormatter();
		CommandLine cmd;

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("utility-name", options);

			System.exit(1);
			return;
		}

		String configFilePath = cmd.getOptionValue("config");
		Boolean isBuilder     = cmd.hasOption("builder");

		
		if(isBuilder) {
			//  *******************
			// ** BUILDING CONFIG ** 
			//  *******************

			BuildingFunction bf = new BuildingFunction(configFilePath, generalSchemaFile, builderSchemaFile, configSchemaFile);
			System.out.println(bf.config.toString(2));
			
			System.exit(0);
			
		} else {
			
			//  ********************
			// ** RUNS CONFIG FILE ** 
			//  ********************
			
	        try {
	        	
	            System.out.println("> Reading configuration JSON object and validating general schema.");
	            config = new JSONModel(configFilePath, generalSchemaFile);
	            	            
	            System.out.println("> Validating configuration schema.");
	            JSONModel.validateConfig(generalSchemaFile, configSchemaFile, config.oModel);
	            
	    		//Adds model fitting constraints if they do not already exist
	    		JSONModel.addEmptyModelFittingConstraints(generalSchemaFile, config.oModel);
	        	
	    		System.out.println("> Builds likelihood object.");
	            CombinedLikelihood l   = new CombinedLikelihood(config.oModel);
	            l.version = version;
	            
	            for(int iComp=0; iComp<l.fittingOrder.size(); iComp++ ) {
	    			l.componentList.get(iComp).setComponentInclusion(true);
	    			l.componentList.get(iComp).setComponentFiting(true);
	    		}
	    		l.updatePacking();
	            
	            System.out.println("> Builds optimizer.");
	            LikelihoodOptimizer lo = new LikelihoodOptimizer(config.oModel);
	            		
	            System.out.println("> Starting optimization.");
	            lo.optimizeLikelihood(l);
	            System.out.println("> Optimization done.");
	            
	            System.exit(0);
	            
	        } catch (java.lang.Exception e) {
	        	e.printStackTrace();
	        	System.exit(1);
	        }
		}
	}
	
	
	public static void testSequenceScorer() {
		
		
		String externalOrder = "AMNTat";
		String alphabetDef = "A-T,M-M,N-N,a-t";
		String seq = "AMNTat";
		
		
		
		
		System.out.println("Alphabet definition         = "+alphabetDef);
		LongSequence.SequenceClass sc = new  LongSequence.SequenceClass(alphabetDef);
		System.out.println("ALPHABET_LOOKUP             = "+sc.getALPHABET_LOOKUP());
		System.out.println("ALPHABET_LOOKUP_REVERSE     = "+sc.getALPHABET_LOOKUP_REVERSE());
		System.out.println("ALPHABET_COMPLEMENT         = "+Misc.formatVector_l(sc.getALPHABET_COMPLEMENT()));
		int alphabetSize = (int) sc.getALPHABET_SIZE();
		System.out.println("ALPHABET_SIZE               = "+alphabetSize);
		int bBits = sc.getbBits();
		System.out.println("bBits (bits/letter)         = "+bBits);
		System.out.println("Bases per long              = "+sc.getBasesPerLong());
		System.out.println("Test sequence               = "+seq);
		System.out.println("Features: ");
		Map<String, List<LongSequence>> features = sc.getFeatures(seq, 1, null);
		for(String k : features.keySet()) {
			System.out.print("  "+k+":  ");
			for(LongSequence l : features.get(k)) {
				System.out.print(l.toString()+", ");
			}
			System.out.println("");
		}
		
		LongSequence ls = sc.build(seq);
		System.out.println("String Complement           = "+ls.getComplement());
		System.out.println("String Reverse              = "+ls.getReverse());
		System.out.println("String getReverseComplement = "+ls.getReverseComplement());
		System.out.println("String Length               = "+ls.getLength());
		List<long[]> swc = ls.getSlidingWindowCodes(2, 0);
		System.out.print("SlidingWindowCodes (L=2)     = [");
		for(long[] l : swc) 
			System.out.print(Misc.formatVector_l(l,",","{","}")+", ");
		System.out.println("]");
		
		ArrayList<ArrayList<double[]>> aal = new ArrayList<ArrayList<double[]>>();
		aal.add(new ArrayList<double[]>());
		aal.add(new ArrayList<double[]>());
		int L=3;
		System.out.println("L                            = "+L);
		
		//Mono betas: (alphabet size) * (window length)
		aal.get(0).add(new double[(int)alphabetSize*L]);
		//Adjacent Di betas: (alphabet size)^2 * (window length -1)
		aal.get(1).add(new double[(L-1)*(int) Math.pow(alphabetSize, 2)]);
		SlidingWindow sw = new SlidingWindow("", "", sc, aal, externalOrder);
		System.out.println("aal[0][0].lenght             = "+aal.get(0).get(0).length);
		System.out.println("L * |alphabet|               = "+(L * alphabetSize));
		System.out.println("aal[1][0].lenght             = "+aal.get(1).get(0).length);
		System.out.println("(L-1) * |alphabet|^2         = "+((L-1)*(int) Math.pow(alphabetSize, 2)));
		int monoCell = (int) Math.pow(2,bBits);
		int diCell = (int) Math.pow(2,2*bBits);
		System.out.println("monoCell                     = "+monoCell);
		System.out.println("diCell                       = "+diCell);
		System.out.println("L*monoCell + (L-1)*diCell    = "+(L*monoCell + (L-1)*diCell));
		System.out.println("Betas                        = "+Misc.formatVector_d(sw.getBeta()));
		System.out.println("Betas.length                 = "+sw.getBeta().length);
		
		ArrayList<ArrayList<Double>> alphas              = sw.slidePN(ls.getValue(), ls.getLength(), 2);
		System.out.println("|alphas|                     = "+alphas.size());
		for(int i=0; i<alphas.size(); i++) 
			System.out.println("   |alphas["+i+"]|               = "+alphas.get(i).size());

		double[] gradient = new double[sw.getBeta().length];
		sw.swGradient(ls.getValue(), ls.getLength(), 1, alphas, gradient);
		System.out.println("|gradient|                   = "+gradient.length);
		int iCurr = 0;
		for(int x=0; x<L; x++) {
			System.out.println("  gradient[mono]["+x+"]          = "+Misc.formatVector_d(Arrays.copyOfRange(gradient, iCurr, iCurr+monoCell), ", ", "{", "}",0).replace("0", " "));
			iCurr += monoCell;
		}
		for(int d=0; d<aal.get(1).size(); d++) {
			for(int x=0; x<L-d-1; x++) {

				System.out.println("  gradient[di]["+d+"]["+x+"]         ={"+Misc.formatVector_d(Arrays.copyOfRange(gradient, iCurr, iCurr+monoCell), ", ", "{", "}",0).replace("0", " ")+",");
				iCurr += monoCell;
				for(int iMono=1; iMono<monoCell; iMono++) {
					System.out.println("                               "+Misc.formatVector_d(Arrays.copyOfRange(gradient, iCurr, iCurr+monoCell), ", ", "{", "}",0).replace("0", " ")+(iMono<monoCell-1?",":"}"));
					iCurr += monoCell;
				}
			}
		}
		System.out.println("Sum(gradient[mono])          = "+(int)Array.sum(Arrays.copyOfRange(gradient, 0, L*monoCell)));
		System.out.println("2 * (|seq|-L+1)*L            = "+2*(seq.length()-L+1)*L);
		System.out.println("Sum(gradient[di])            = "+(int)Array.sum(Arrays.copyOfRange(gradient,  L*monoCell, gradient.length)));
		System.out.println("2 * (|seq|-L+1)*(L-1)        = "+2*(seq.length()-L+1)*(L-1));
		
		double[] convGradient = sw.convertGradient(gradient, externalOrder);
		System.out.println();
		System.out.println("External Ordering            = "+externalOrder);
		System.out.println("|converted gradient|         = "+convGradient.length);
		iCurr = 0;
		for(int x=0; x<L; x++) {
			System.out.println(" converted gradient[mono]["+x+"] = "+Misc.formatVector_d(Arrays.copyOfRange(convGradient, iCurr, iCurr+alphabetSize), ", ", "{", "}",0).replace("0", " "));
			iCurr += alphabetSize;
		}
		for(int d=0; d<aal.get(1).size(); d++) {
			for(int x=0; x<L-d-1; x++) {
				System.out.println(" converted gradient[di]["+d+"]["+x+"] ={"+Misc.formatVector_d(Arrays.copyOfRange(convGradient, iCurr, iCurr+alphabetSize), ", ", "{", "}",0).replace("0", " ")+",");
				iCurr += alphabetSize;
				for(int iMono=1; iMono<alphabetSize; iMono++) {
					System.out.println("                                "+Misc.formatVector_d(Arrays.copyOfRange(convGradient, iCurr, iCurr+alphabetSize), ", ", "{", "}",0).replace("0", " ")+(iMono<monoCell-1?",":"}"));
					iCurr += alphabetSize;
				}
			}
		}

		System.out.println("Sum(converted gradient[mono])= "+(int)Array.sum(Arrays.copyOfRange(convGradient, 0, L*alphabetSize)));
		System.out.println("Sum(converted gradient[di])  = "+(int)Array.sum(Arrays.copyOfRange(convGradient,  L*alphabetSize, gradient.length)));

		
	}

	
	
}
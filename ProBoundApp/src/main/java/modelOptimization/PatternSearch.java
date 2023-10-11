package modelOptimization;

import java.util.Random;

import org.json.JSONObject;

import base.Array;

public class PatternSearch extends Optimizer {

	private boolean isVerbose, useRandom;
	private int totFeatures;
	private double d0, theta, epsilon;

	
	public PatternSearch(JSONObject settings, /*double d0, double dTol, double theta, boolean useRandom,*/ 
			boolean isVerbose) {
		this.d0					= settings.getDouble("initStep");
		this.epsilon			= settings.getDouble("convergence");
		this.theta				= settings.getDouble("theta");
		this.useRandom			= settings.getBoolean("randomAxis");
		this.isVerbose			= isVerbose;
	}
	
	@Override
	public JSONObject saveToJSON_settings() {
		JSONObject out = new JSONObject();
		out.put("initStep",    d0);
		out.put("convergence", epsilon);
		out.put("theta",       theta);
		out.put("randomAxis",  useRandom);
		return out;
		
	}
	

	
	@Override
	public boolean optimize(LossFunction lf) throws Exception{
		
		//this.lf                 = lf;
		totFeatures             = lf.nParameters;
		boolean hasImproved		= false;						//Has the current step improved the estimate?
		int totalSteps			= 0;							//total number of steps taken
		int stepPos;
		int functionCalls		= 1;
		double d 				= d0;							//step size
		double delta			= 0;							//distance moved per step
		double forwardStep		= 0;
		double reverseStep		= 0;
		double tStart			= System.currentTimeMillis();
		double[] prevLoc;
		double[] currBetas		= null;
		double[] forwardBetas	= null;
		double[] reverseBetas	= null;
		double currBestFVal 	= 0;
		int[] featureOrder		= null;
		Random generator		= new Random();
		
		currBetas = Array.clone(lf.parameters);
		lf.lossFunction_updateValue();
		currBestFVal	= lf.value;
		if (isVerbose){
			System.out.println("Starting Function Value: "+currBestFVal);
			System.out.println("Iterations   Fnc. Calls           Likelihood"
					+ "       Distance Moved            Step Size");
		}
		prevLoc			= Array.clone(currBetas);
		
		while (true) {
			if (d < epsilon) {											//check to see if convergence has been reached
				tStart		= (System.currentTimeMillis()-tStart)/1000;
				lf.lossFunction_reportPosition();
				return true;
			} else {												//if not, begin mutation steps
				hasImproved = false;								//reset	
				//Loop over all nucleotide features
				if (useRandom)	featureOrder = getRandomPermutation(totFeatures, generator);
				for (int phi=0; phi<totFeatures; phi++) {
					stepPos	= (useRandom) ? featureOrder[phi] : phi;
					//Check to see if position has already been dealt with

					//Compute forward and reverse betas (and symmetrize)
					forwardBetas     = Array.clone(currBetas);
					forwardBetas[stepPos] += d;
					lf.lossFunction_setParameters(forwardBetas);
					lf.lossFunction_updateValue();
					forwardStep	     = lf.value;
					
					reverseBetas     = Array.clone(currBetas);
					reverseBetas[stepPos] -= d;
					lf.lossFunction_setParameters(reverseBetas);
					lf.lossFunction_updateValue();
					reverseStep	     = lf.value;
					functionCalls+= 2;
					
					
					//Forward direction improved log likelihood
					if((currBestFVal-forwardStep)/currBestFVal > 5e-16) {
						//Reverse direction is even better
						if ((forwardStep-reverseStep)/forwardStep > 5e-16) {
							currBetas	= reverseBetas;
							currBestFVal= reverseStep;
							hasImproved = true;
						} else {									//forward direction is better
							currBetas	= forwardBetas;
							currBestFVal= forwardStep;
							hasImproved	= true;
						}
					} else if ((currBestFVal-reverseStep)/currBestFVal > 5e-16) { //reverse direction improved log likelihood
						currBetas	= reverseBetas;
						currBestFVal= reverseStep;
						hasImproved	= true;
					}
					
				}				
				if (hasImproved) {									//has there been an reduction in the function value?
					totalSteps++;
					delta 		= Array.dist(prevLoc, currBetas);
					prevLoc		= Array.clone(currBetas);
					if (isVerbose) {
						printStep(totalSteps, functionCalls, currBestFVal, delta, d);
					}
					lf.lossFunction_reportPosition();
				} else {
					d = d*theta;									//and reduce step size
				}
			}
		}	
	}
	
	public static int[] getRandomPermutation (int length, Random generator){
	    int idx		= 0;
	    int swap	= 0;
	    int[] array = new int[length];
	 
	    // initialize array and fill it with {0,1,2...}
	    for(int i = 0; i < array.length; i++) {
	    	array[i] = i;
	    }

	    for(int i = 0; i<length; i++){
	        // randomly chosen position in array whose element
	        // will be swapped with the element in position i
	        // note that when i = 0, any position can chosen (0 thru length-1)
	        // when i = 1, only positions 1 through length -1
	        idx = i+generator.nextInt(length-i);

	        // perform swap
	        swap		= array[i];
	        array[i] 	= array[idx];
	        array[idx] 	= swap;
	    }                       
	    return array;
	}

}

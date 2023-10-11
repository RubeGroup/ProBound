package modelOptimization;

import proBoundApp.Array;
import proBoundApp.Misc;

abstract public class LossFunction {
	
	public int nParameters;      //Number of parameters
	public double[] parameters;  // Current value of parameter. Set by setParameters(double[]);
	public double value;         // Current value of the loss function. Computed by updateValue()
	public double[] gradient;    // Current value of the gradient of the loss function. Computed by updateGradient();
	public boolean isNan;        // Indicates that the evaluation has failed. Indicates that an evaluation has failed.
	public double[] lastReported;// Contains the the last reported parameters.
	
	protected boolean computeVariance = false; //Indicates if the variance of the value and gradient should be computed
	public double valueVariance;               //Variance of the function value
	public double[] gradientVariance;          //Variance of the components in the gradient vector
		
	abstract public void lossFunction_setParameters(double[] in);  // Sets the parameters
	abstract public void lossFunction_updateValue();               // Computes the value of the loss function
	abstract public void lossFunction_updateGradient();            // Computes the gradient
	abstract public void lossFunction_reportPosition();            // Is called after the optimizer takes as step. Used for logging
	abstract public void lossFunction_nextBatch(int nData);        // Gets a new batch of nData data points to use for fuction/gradient evaluation.
	abstract public long lossFunction_getDataSize();	           // Returns the number of datapoints.
	abstract public int  lossFunction_getBatchSize();	           // Returns the number of datapoints in the batch.
	
	// Sets whether or not the variance of the function value should be computed.
	public void lossFunction_setComputeVariance(boolean use) {
		
		computeVariance = use;
		
		if(use) {
			valueVariance    = 0;
			gradientVariance = new double[nParameters];
		} else {
			valueVariance    = 0;
			gradientVariance = null;
		}
	}
	
	
	public double[] gradientFiniteDifferences(double stepSize) {
		int totFeatures		= parameters.length;
		double[] startLoc   = Array.clone(parameters) ;
		double forwardDifference, reverseDifference;
		double[] modBetas;
		double[] fdGradient	= new double[totFeatures];
		
		//compute function values
		for (int i=0; i<totFeatures; i++) {
			try {
				modBetas			 = Array.clone(startLoc);
				modBetas[i]			+= stepSize;
				lossFunction_setParameters(modBetas);
				lossFunction_updateValue();
				forwardDifference    = value;
				modBetas[i]         -= 2*stepSize;
				lossFunction_setParameters(modBetas);
				lossFunction_updateValue();
				reverseDifference    = value;
				fdGradient[i]        = (forwardDifference-reverseDifference)/(2*stepSize);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		lossFunction_setParameters(startLoc);

		return fdGradient;
	}
	
	public double[][] hessianFiniteDifferences(double stepSize) {
		double[] startLoc = Array.clone(parameters);
		int totFeatures 			= parameters.length;
		double[] baseBetas, modBetas;
		double[][] forwardDifference= new double[totFeatures][totFeatures];
		double[][] reverseDifference= new double[totFeatures][totFeatures];
		double[][] fdHessian		= new double[totFeatures][totFeatures];
		
		baseBetas = Array.clone(startLoc);
		//compute gradients
		for (int i=0; i<totFeatures; i++) {
			try {
				modBetas			= Array.clone(baseBetas);
				modBetas[i]			+= stepSize;
				lossFunction_setParameters(modBetas);
				lossFunction_updateGradient();
				forwardDifference[i]= gradient;
				modBetas[i]			-= 2*stepSize;
				lossFunction_setParameters(modBetas);
				lossFunction_updateGradient();
				reverseDifference[i]= gradient;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		//Find finite differences (forward/backward FD)
		for (int i=0; i<totFeatures; i++) {
			for (int j=0; j<totFeatures; j++) {
				fdHessian[i][j] = (forwardDifference[j][i]-reverseDifference[j][i])/(4*stepSize) + 
							      (forwardDifference[i][j]-reverseDifference[i][j])/(4*stepSize);
			}
		}
		lossFunction_setParameters(startLoc);
		return fdHessian;
	}
	
	public void checkVariance() {

		//Testing so sub-batches that the expected value is constant across batch sizes and that the variance scales like 1/sqrt(n)
		System.out.println("");
		System.out.println("Starting to check variance and batches...");
		System.out.println("=========================================");
		System.out.println("Number of datapoints: "+lossFunction_getDataSize());
		int nDec = (int) Math.ceil(Math.log10(lossFunction_getDataSize()));
		int[] depth = new int[nDec];
		for(int iDec=0; iDec<nDec; iDec++)
			depth[iDec] =(int) Math.pow(10, iDec);

		System.out.println("");
		System.out.println("Checking how the function value (and its variance across probes) changes with batch size.");
		System.out.println("-----------------------------------------------------------------------------------------");
		lossFunction_setComputeVariance(true);
		lossFunction_nextBatch(-1);
		lossFunction_updateValue();
		double fullValue = value;
		System.out.println("Full loss function:                        \t"+String.format("%.4f", fullValue));
		lossFunction_nextBatch((int)lossFunction_getDataSize()-1);
		lossFunction_updateValue();
		System.out.println("Loss function after dropping one datapoint:\t"+String.format("%.4f", value));
		System.out.println("Absolute difference (should be small):     \t"+String.format("%.4e", Math.abs(fullValue-value))+"\n");
		
		int nTest = 100;
		System.out.println("Computing function values using varying batch size. n="+nTest+" replicates for each.\n");
		System.out.println("#Batch  | Mean across batches +/- SD | SD * sqrt(n_batch)  | z-score of mean | sqrt(Mean value variance)");
		System.out.println("--------+----------------------------+---------------------+-----------------+--------------------------");
		for(int iD=0; iD<depth.length; iD++){
			int d=depth[iD];
		
			double[] funcValues = new double[nTest];
			double[] varFuncValues = new double[nTest];
			for(int iTest=0; iTest<nTest; iTest++) {
				lossFunction_nextBatch(d);
				lossFunction_updateValue();
				funcValues[iTest] = value;
				varFuncValues[iTest] = valueVariance;
			}
			
			double funcMean = Array.mean(funcValues);
			double sigma    = Math.sqrt(Array.dotProduct(funcValues,funcValues)/nTest - Math.pow(funcMean, 2));

			System.out.println(
					String.format("%7d",d)+" | "
							+ String.format("%10.4f",funcMean) + " +/- " + String.format("%9.4f",sigma) 
							+ "   | \t\t" + String.format("%10.4f",sigma*Math.sqrt(d)) 
							+ " | \t"+String.format("%9.4f",Math.abs(funcMean-fullValue)/(sigma/Math.sqrt(nTest)))
							+ "    | \t"+String.format("%9.4f",Math.sqrt(Array.mean(varFuncValues))));
		}
		System.out.println("");
		System.out.println("    What to expect:");
		System.out.println("    \"Mean across batches\" should be constant                  => Batch size does not bias function value.");
		System.out.println("    \"SD *sqrt(n_batch)\" should be constant.                   => Stochastic noise decreasese like 1/sqrt(n)");
		System.out.println("    \"z-score of mean\" should follow N(0,1).                   => We understand distribution.");
		System.out.println("    \"sqrt(Men value variance)\"  sould match SD * sqrt(nBatch) => in-batch across-batch variance consistent.\n");

		
		
		System.out.println("Checking how the gradient value (and its variance across probes) changes with batch size.");
		System.out.println("-----------------------------------------------------------------------------------------");
		lossFunction_nextBatch(-1);
		lossFunction_updateGradient();
		double[] fullGrad = Array.clone(gradient);
		System.out.println("Full gradient:                                  "+Misc.formatVector_d(fullGrad,",","{","}",6));
		lossFunction_nextBatch((int)lossFunction_getDataSize()-1);
		lossFunction_updateGradient();
		System.out.println("Gradient function after dropping one datapoint: "+Misc.formatVector_d(gradient,",","{","}",6));
		System.out.println("Mean(|full-drop1|):                             "+String.format("%.4e", Array.mean(Array.abs(Array.subtract(fullGrad, gradient)))));
		
		System.out.println("");
		System.out.println("#Batch  | Mean |grad| over batches +/- SD | SD * sqrt(n_batch) | z-score of mean | sqrt(Mean gradient variance)");
		System.out.println("--------+---------------------------------+--------------------+-----------------+-----------------------------");

		double[] gradSum    = null;
		double[] gradSum2   = null;
		double[] grandSumSE = null;
		double[] varGradSum = null;
		for(int d : depth) {
			
			//Computes the mean gradient, the component-wise mean of gradient, gradient^2, (gradient - E[Gradient])^2)
			gradSum    = null;
			gradSum2   = null;
			grandSumSE = null;
			varGradSum = null;
			for(int iTest=0; iTest<nTest; iTest++) {
				lossFunction_nextBatch(d);
				lossFunction_updateGradient();
				double[] tempGrad    = gradient;
				double[] tempGrad2   = Array.multiply(gradient, gradient);
				double[] tempSE      = Array.pow(Array.subtract(gradient, fullGrad), 2);
				double[] tempVarGrad = gradientVariance;
				gradSum    = gradSum   ==null ? tempGrad    : Array.add(gradSum,    tempGrad);
				gradSum2   = gradSum2  ==null ? tempGrad2   : Array.add(gradSum2,   tempGrad2);
				grandSumSE = grandSumSE==null ? tempSE      : Array.add(grandSumSE, tempSE);
				varGradSum = varGradSum==null ? tempVarGrad : Array.add(varGradSum, tempVarGrad);
			}
			double[] gradMean    = Array.scalarMultiply(gradSum,    1./nTest);
			double[] grad2Mean   = Array.scalarMultiply(gradSum2,   1./nTest);
			double[] varGradMean = Array.scalarMultiply(varGradSum, 1./nTest);
			double[] gradSigma   = Array.sqrt(Array.subtract(grad2Mean, Array.multiply(gradMean, gradMean)));
			//double[] gradMSE     = Array.scalarMultiply(grandSumSE, 1./nTest);
			double[] absDiff     = Array.abs(Array.subtract(fullGrad, gradMean));
			
			System.out.println(
					String.format("%7d",d)+" | "
							+ String.format("%10.4f",Array.mean(Array.abs(gradMean))) + " +/- " + String.format("%9.4f",Array.mean(gradSigma)) 
							+ "        |         " + String.format("%10.4f",Array.mean(gradSigma)*Math.sqrt(d)) 
							+ " |    "+String.format("%9.4f",Array.mean(Array.scalarMultiply(Array.divide(absDiff, gradSigma), Math.sqrt(nTest))))
							+ "    | \t"+String.format("%9.4f",Array.mean(Array.sqrt(varGradMean))));
			
		}
	}
}

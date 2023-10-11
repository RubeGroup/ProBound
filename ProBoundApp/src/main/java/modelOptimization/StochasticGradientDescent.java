package modelOptimization;

import org.json.JSONObject;

import modelComponents.ModelComponent;

public class StochasticGradientDescent extends Optimizer {

	String method;
	int batchsize, maxIters;
	double learningRate, convergence;
	double[] adaptive;
	boolean isVerbose;
	
	public StochasticGradientDescent(JSONObject settings, boolean verbose) {
		
		isVerbose = verbose;
		
		this.method       = settings.getString("method");
		this.batchsize    = settings.getInt("batchSize");
		this.learningRate = settings.getDouble("learningRate");
		this.convergence  = settings.getDouble("convergence");
		this.maxIters     = settings.getInt("maxIters");
		this.adaptive     =  ModelComponent.readFromJSON_d(settings.getJSONArray("adaptive"));
		
		System.out.println("StochasticGradientDescent object created.");
	}
	
	@Override
	public JSONObject saveToJSON_settings() {
		
		JSONObject out = new JSONObject();
		out.put("method",       method);
		out.put("batchSize",    batchsize);
		out.put("learningRate", learningRate);
		out.put("convergence",  convergence);
		out.put("maxIters",     maxIters);
		out.put("adaptive",     ModelComponent.JSONArrayConvert_d(adaptive));
		return out;
		
	}
		
	@Override
	public boolean optimize(LossFunction l) throws Exception {
		return true;
	}

}

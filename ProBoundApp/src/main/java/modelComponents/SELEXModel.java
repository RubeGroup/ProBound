package modelComponents;

import java.util.ArrayList;
import java.util.Hashtable;

import org.json.JSONObject;

public class SELEXModel extends EnrichmentModel{
	
	public boolean bindingSaturation;
	public boolean trySaturation;
	//double[] rounds;
	SELEXModel(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModes, 
			ArrayList<BindingModeInteraction> allInteractions) {
		
		super(config, iExpIn, inTable, allBindingModes, allInteractions);
		
		componentName = "SELEX enrichment model "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");
		
		bindingSaturation = false;
		trySaturation     = false;
		
		readFromJSON_settings(config);
		readFromJSON_constraints(config);

		allocateParameters();
		
		setFreezingLevel(0);
		
		variationName       = null;
		variationsOptimized = false;

	}
	
	public void allocateParameters() {
		
		return;
		
	}
	
	
	@Override
	public void activationAdjustment(double weight) {
		
		standardActivityAdjustment(weight);
		
	}
	
	@Override
	public Hashtable<String,JSONObject> getVariations(JSONObject currentModel, ArrayList<ModelComponent> componentList) {
		if(variationsOptimized)
			return null;
		
		//Initial model shift
		/////////////////
		Hashtable<String,JSONObject> newVariations = new Hashtable<String,JSONObject>();
		String testLabel = "saturation";
		if(!variationDescription.containsKey(testLabel)) {
			variationDescription.put(testLabel, "Model with saturation.");
			JSONObject modelVariation = clone_JSON_O(currentModel);
			modelVariation.getJSONObject("modelSettings").getJSONArray("enrichmentModel").getJSONObject(iComp).put("bindingSaturation", true);
			newVariations.put(testLabel, modelVariation);
		}
		
		variationsOptimized = true;
		
		if(newVariations.size()>0)
			return newVariations;
		
		//No variation possible.
		return null;

	}
	
	@Override
	public void seed_component(JSONObject config) {

		return;
	}
	
	@Override
	int packModel_component(JSONObject packing, int iFirst) {
		
		String coefficientKey = "packing";
		addEmptyJSON_component_O(packing, coefficientKey, componentKey, iComp);
		int iCurr = iFirst;
		return iCurr;
	}
	

	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		
		return;
	}
	
	@Override
	void saveToJSON_settings(JSONObject out) {
		
		super.saveToJSON_settings(out);
		
		String coefficientKey = "modelSettings";
		JSONObject oEnr   =                  out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oEnr.put("bindingSaturation",  bindingSaturation);

				
	}
	
	@Override
	void readFromJSON_settings(JSONObject in) {
		
		String coefficientKey = "modelSettings";
		super.readFromJSON_settings(in);
		JSONObject oSettEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		bindingSaturation     = oSettEnr.getBoolean("bindingSaturation");

	}
	
	
	@Override
	void saveToJSON_constraints(JSONObject out) {
		
		String coefficientKey = "modelFittingConstraints";
		addEmptyJSON_component_O(out, coefficientKey, componentKey, iComp);
		
		JSONObject oEnr       = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);

		oEnr.put("trySaturation",             trySaturation);
		
	}
	
	@Override
	void readFromJSON_constraints(JSONObject in) {
		
		String coefficientKey = "modelFittingConstraints";
		if(in.has(coefficientKey)) {
			JSONObject oConsEnr   = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
			trySaturation          = oConsEnr.has("trySaturation") ? oConsEnr.getBoolean("trySaturation") : false;		
		}
	}
	
	@Override
	public void regularizationToJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		
	}
	
	@Override
	void saveToJSON_parameters(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		
	}

	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {
		
		//JSONObject oEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iExp);
	
	}
	
	public void updateDeltaKappaRI(double[] deltaKappaRI, double[] alphaRI, int nRounds) {
		for(int r=0; r<nRounds; r++) {					
			deltaKappaRI[r] = alphaRI[r];
			if(bindingSaturation)
				deltaKappaRI[r] /= (1+alphaRI[r]);
		}
	}
	
	public void updateNablaF1(double[] nablaF1, double[] nablaN, double[] alphaRI, int nColumns) {
		for(int r=0; r<nColumns; r++) {
			nablaF1[r]      = nablaN[r] * (  1.0/alphaRI[r]     );
			if(bindingSaturation)
				nablaF1[r] += nablaN[r] * ( -1.0/(1+alphaRI[r]) );		
		}
	}
	
	//Code for saving a contribution to the EnrichmentModel gradient.
	/////////////////////////////////////////////////////////////////
	public enrichmentGradient newGradient() {
		return new SELEXGradient();
	}

	public class SELEXGradient extends enrichmentGradient {

		SELEXGradient() { super(); }

		void addEnrichmentGradientContribution(double[] nablaN, double[] alphaRI) { }

		void saveToJSON(JSONObject outObject, String coefficientKey) { }
		void saveToJSON_squared(JSONObject outObject, String coefficientKey) { }

	}
}

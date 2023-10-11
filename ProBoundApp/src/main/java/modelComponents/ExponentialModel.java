package modelComponents;

import java.util.ArrayList;
import java.util.Hashtable;

import org.json.JSONObject;

public class ExponentialModel extends EnrichmentModel {

	ExponentialModel(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModes, 
			ArrayList<BindingModeInteraction> allInteractions) {
		
		super(config, iExpIn, inTable, allBindingModes, allInteractions);
		
		iComp          = iExpIn;
		componentName = "Exponential enrichment model "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");

		modelType = "Exponential";
		
		readFromJSON_settings(config);
		readFromJSON_constraints(config);

		allocateParameters();
		
		setFreezingLevel(0);
		
		variationName       = null;
		variationsOptimized = true;

	}
	
	public void allocateParameters() {
		
		return;
		
	}
	
	
	@Override
	public void activationAdjustment(double weight) {
		
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
 
	}
	
	@Override
	void readFromJSON_settings(JSONObject in) {

		super.readFromJSON_settings(in);

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
	}
	
	@Override
	void readFromJSON_constraints(JSONObject in) {
		
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
		
	}

	public void updateDeltaKappaRI(double[] deltaKappaRI, double[] alphaRI, int nRounds) {
		for(int r=0; r<nRounds; r++) 			
			deltaKappaRI[r] = Math.exp(alphaRI[r]);
	}
	
	public void updateNablaF1(double[] nablaF1, double[] nablaN, double[] alphaRI, int nColumns) {
		for(int r=0; r<nColumns; r++)
			nablaF1[r] = nablaN[r];
	}
	
	//Code for saving a contribution to the EnrichmentModel gradient.
	/////////////////////////////////////////////////////////////////
	public enrichmentGradient newGradient() {
		return new ExponentialGradient();
	}

	public class ExponentialGradient extends enrichmentGradient {

		ExponentialGradient() { super(); }

		void addEnrichmentGradientContribution(double[] nablaN, double[] alphaRI) { }

		void saveToJSON(JSONObject outObject, String coefficientKey) { }
		void saveToJSON_squared(JSONObject outObject, String coefficientKey) { }

	}
}

package modelComponents;

import java.util.ArrayList;
import java.util.Hashtable;

import org.json.JSONObject;

import proBoundApp.Misc;


public class InputFreeFraction extends EnrichmentModel {
	
	// Enrichment model describing a multi-band experiment with columns 
	// Col 0 = Input                 = 1 + Z_2 + Z_3 + ...
	// Col 1 = Free                  = 1
	// Col 2 = Fraction 1 (Monomer)  = Z_2
	// Col 3 = Fraction 2 (Dimer)    = Z_3
	// ...
	
	InputFreeFraction(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModes, 
			ArrayList<BindingModeInteraction> allInteractions) {
		
		super(config, iExpIn, inTable, allBindingModes, allInteractions);
		
		componentName = "Input-Free-Fraction enrichment model "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");

		modelType = "RhoGamma";
		
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

		return iFirst;
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
		addEmptyJSON_component_O(out, coefficientKey, componentKey, iComp);
		
	}
	
	@Override
	void readFromJSON_constraints(JSONObject in) {
	
		return;
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
		
		return;
	
	}

	public void updateDeltaKappaRI(double[] deltaKappaRI, double[] alphaRI, int nRounds) {
		
		// Z_tot = 1 + Z_2 + Z_3
		// Col 0 = Input                 = 1
		// Col 1 = Free                  = 1   / Z_tot
		// Col 2 = Fraction 1 (Monomer)  = Z_2 / Z_tot
		// Col 3 = Fraction 2 (Dimer)    = Z_3 / Z_tot
		
		double Ztot=1.0;
		for(int r=2; r<nColumns; r++)
			Ztot += alphaRI[r];

		
		
		deltaKappaRI[0] = 1;
		deltaKappaRI[1] = 1/Ztot;
		for(int r=2; r<nColumns; r++) 
			deltaKappaRI[r]  = alphaRI[r]/Ztot;

	}
	
	public void updateNablaF1(double[] nablaF1, double[] nablaN, double[] alphaRI, int nRounds) {
		
		// Col 0 = Input                 dK = 0
		// Col 1 = Free                  dK = -nablaNSum / Ztot
		// Col 2 = Fraction 1 (Monomer)  dK = -nablaNSum / Ztot + nablaN2 / alpha2
		// Col 3 = Fraction 2 (Dimer)    dK = -nablaNSum / Ztot + nablaN2 / alpha2

		//Computes 1+Z_2+Z_3+...
		double Ztot = 0;
		double nablaNSum = 0.0;

		for(int r=1; r<nColumns; r++) {
			Ztot      += r==1 ? 1 : alphaRI[r];
			nablaNSum += nablaN[r];
		}
		
		
		//Computes d(log L)/d(alpha_r) = sum_t delta n_t d(kappa_t)/d(alpha_t) / kappa_t
		nablaF1[0] = 0;
		nablaF1[1] = 0;
		for(int r=2; r<nColumns; r++)
			nablaF1[r] = nablaN[r]/alphaRI[r] - nablaNSum/Ztot;
		
	}
	
	//Code for saving a contribution to the EnrichmentModel gradient.
	/////////////////////////////////////////////////////////////////
	public enrichmentGradient newGradient() {
		return new imputFreeFractionGradient();
		
	}
	
	public class imputFreeFractionGradient extends enrichmentGradient {
		
		imputFreeFractionGradient() {
			super();
		}
		
		void addEnrichmentGradientContribution(double[] nablaN, double[] alphaRI) {
			return;
			
		}
		
		void saveToJSON(JSONObject outObject, String coefficientKey) {
			return;
		}
		
		void saveToJSON_squared(JSONObject outObject, String coefficientKey) {
			return;
		}
		
	}
	
	static public void printJSONObjectCoefficients(JSONObject model, String coefficientKey, int iExp) {

		return;
	}
}

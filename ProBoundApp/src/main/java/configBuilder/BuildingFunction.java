package configBuilder;

import org.json.*;

import modelComponents.ModelComponent;
import modelOptimization.CombinedLikelihood;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.sql.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import proBoundApp.JSONModel;

public class BuildingFunction {

	// COMMAND OVERVIEW:
	// addTable(countTableFile, variableRegionLength, nColumns, ...)   //Explicitly adds a single count table. 
	// addTableDB(count_table_id, ...)                                 //Adds a single count table from CELLX
	// addTableDBs(count_table_ids)	                                   //Adds multiple tables forom CellX
	// addSELEXTableDBs(count_table_ids, ...)                          //Adds multiple count tables.
	// addSELEX(...)                                                   //Explicitly adds a SELEX-experiment enrichment model
	// addSELEXDB(count_table_id, ...)                                 //Adds a single SELEX experiment using the CellX database.
	// addSELEXTableDB(count_table_id, ...)                            //Adds both a single count table and a single SELEX experiment.
	// addBindingMode(size, ...)                                       //Adds a binding mode
	// addNS()                                                         //Adds a binding mode with NS binding mode
	// addInteraction(bindingModes, ...)                               //Adds a binding mode interactions
	// output(outputPath, baseName, ...)                               //Adds output.
	// outputDB(fit_id)                                                //Adds output using CellX
	// optimizerSetting(...)                                           //Adds optimizer settings.
	// lbfgsSettings(...)                                              //Adds LBFGS settings
	// slbfgsSettings(...)                                             //Adds S-LBFGS settings.
	// patternSearchSettings(...)                                      //Adds pattern search settings
	// basicConstraints(...)                                           //Adds basic model fitting constraints (i.e. addBindingModesSequentially).
	// enrichmentModelConstraints(...)                                 //Adds enrichment model constraints.
	// bindingModeConstraints(...)                                     //Adds binding mode constraints.
	// interactionConstraints(...)                                     //Adds interactions constraints.
	// NRLBConstraints(...)                                            //Adds default NRLB constraints.
	// symmetry(...)                                                   //Sets the symmetry of a binding mode.
	// bindingModeSeed(...)                                            //Adds binding mode seed.
	// enrichmentModelSeed(...)                                        //Adds enrichment model seed.
	// setAlphabet(...)                                                //Sets the alphabet
	
	// USAGE FOR NRLB:
	// IDEA:
	// new_fitEntry.py -r ProBound -t 2 -v size=10 -c standard.nsDi_NRLB.json

	private String cellx_dir;
	public JSONObject config;
	Map<String, String> cellx_config;

	public BuildingFunction() {


		//Reads CELLX environment.
		if(System.getenv("CELLX_DIR")==null) {
			cellx_dir    = null;
			cellx_config = null;
		} else {
			cellx_dir   = System.getenv("CELLX_DIR");
			String conf = cellx_dir+"/code/config.sh";
			
			if(!(new File(conf)).exists())
				throw new java.lang.RuntimeException("ERROR: The CELLX config file "+conf+" does not exist.");
			try {
				Process p = new ProcessBuilder("/bin/bash", "-c", ". "+conf+"; env | grep ^CELLX_").start();
				p.waitFor();
				BufferedReader reader = 
						new BufferedReader(new InputStreamReader(p.getInputStream()));

				cellx_config = new HashMap<String, String>();
				String line = "";           
				while ((line = reader.readLine())!= null) {
					String[] d = line.split("=");
					if(d.length == 2)
						cellx_config.put(d[0], d[1]);
				}

			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}

	}
	
	public BuildingFunction(JSONArray buildSteps, String generalSchemaFile, String builderSchemaFile, String configSchemaFile) {
		this();
		buildConfig(buildSteps, generalSchemaFile, builderSchemaFile, configSchemaFile);
	}
	
	public BuildingFunction(String jsonPath, String generalSchemaFile, String builderSchemaFile, String configSchemaFile) {
		
		this();
		
		JSONArray parsedJSONArray = null;

		try {

			File jsonFile                    = new File(jsonPath);
			InputStream is                   = new FileInputStream(jsonFile);
			BufferedReader streamReader      = new BufferedReader(new InputStreamReader(is, "UTF-8"));
			StringBuilder responseStrBuilder = new StringBuilder();

			String inputStr;
			while ((inputStr = streamReader.readLine()) != null)
				responseStrBuilder.append(inputStr);

			parsedJSONArray = new JSONArray(responseStrBuilder.toString());

			streamReader.close();

		} catch (JSONException e) {
			throw new java.lang.RuntimeException("JSONException reading input JSON file '"+jsonPath+"': "+e.getMessage(), e);

		} catch (java.io.FileNotFoundException e) {
			
			System.err.println("ERROR: Can't open file:");
			System.err.println(e.getMessage());
			System.exit(1);
			
		} catch (java.lang.Exception e) {
			throw new java.lang.RuntimeException("Problem opening and reading the file '"+jsonPath+"'.", e);
		}

		buildConfig(parsedJSONArray, generalSchemaFile, builderSchemaFile, configSchemaFile);
	}
	
	public void buildConfig(JSONArray builtSteps, String generalSchemaFile, String builderSchemaFile, String configSchemaFile) {
		
		JSONModel.validateSchemaFile_A(builderSchemaFile, builtSteps);
		
		config = new JSONObject();
		
		for(int i=0; i<builtSteps.length(); i++) {
			JSONObject step = builtSteps.getJSONObject(i);
			String f = step.getString("function");
			if(f.equals("addTable"))
				addTable(step);
//			else if(f.equals("addTableDB"))
//				addTableDB(step);
//			else if(f.equals("addTableDBs"))
//				addTableDBs(step);
			else if(f.equals("addSELEX"))
				addSELEX(step);
//			else if(f.equals("addSELEXDB"))
//				addSELEXDB(step);
//			else if(f.equals("addSELEXTableDB"))
//				addSELEXTableDB(step);
//			else if(f.equals("addSELEXTableDBs"))
//				addSELEXTableDBs(step);
			else if(f.equals("addBindingMode"))
				addBindingMode(step);
			else if(f.equals("addNS"))
				addNS(step);
			else if(f.equals("addInteraction"))
				addInteraction(step);
//			else if(f.equals("outputDB"))
//				outputDB(step);
			else if(f.equals("output"))
				output(step);
			else if(f.equals("optimizerSetting"))
				optimizerSetting(step);
			else if(f.equals("lbfgsSettings"))
				lbfgsSettings(step);
			else if(f.equals("patternSearchSettings"))
				patternSearchSettings(step);
			else if(f.equals("slbfgsSettings"))
				slbfgsSettings(step);
			else if(f.equals("basicConstraints"))
				basicConstraints(step);
			else if(f.equals("enrichmentModelConstraints"))
				enrichmentModelConstraints(step);
			else if(f.equals("bindingModeConstraints"))
				bindingModeConstraints(step);
			else if(f.equals("interactionConstraints"))
				interactionConstraints(step);
			else if(f.equals("NRLBConstraints"))
				NRLBConstraints(step);
			else if(f.equals("symmetry"))
				symmetry(step);
			else if(f.equals("bindingModeSeed"))
				bindingModeSeed(step);
			else if(f.equals("enrichmentModelSeed"))
				enrichmentModelSeed(step);
			else if(f.equals("setAlphabet"))
				setAlphabet(step);
			else
				throw new java.lang.RuntimeException("ERROR: Invalid builer function: "+step.getString("function"));
		}
		
		JSONModel.validateSchemaFile_O(generalSchemaFile, config);
		JSONModel.validateSchemaFile_O(configSchemaFile, config);

	}
	
	private Connection getDBConnection() {
		Connection c = null;

		if(cellx_dir==null)
			throw new java.lang.RuntimeException("Cannot connect to the database since CELLX_DIR has not been set.");
	
		try {
			Class.forName("org.postgresql.Driver");
			c = DriverManager
					.getConnection("jdbc:postgresql://"+cellx_config.get("CELLX_SQL_IP")+":"+cellx_config.get("CELLX_SQL_PORT")+"/"+cellx_config.get("CELLX_SQL_DATABASE"),
							cellx_config.get("CELLX_SQL_USER"), cellx_config.get("CELLX_SQL_PASSWORD"));

		} catch (Exception e) {
			e.printStackTrace();
			System.err.println(e.getClass().getName()+": "+e.getMessage());
			System.exit(1);
		}

		return c;
	}

	JSONObject createAndReturnEntry_O(String key1, String key2) {
		
		if(!config.has(key1))
			config.put(key1, new JSONObject());
		JSONObject o1 = config.getJSONObject(key1);
		
		if(!o1.has(key2))
			o1.put(key2, new JSONObject());
		JSONObject o2 = o1.getJSONObject(key2);
		
		return o2;
		
	}
	
	JSONObject createAndReturnEntry_A(String key1, String key2, int iCoeff) {
		
		if(!config.has(key1))
			config.put(key1, new JSONObject());
		JSONObject o1 = config.getJSONObject(key1);
		
		if(!o1.has(key2))
			o1.put(key2, new JSONArray());
		JSONArray a2 = o1.getJSONArray(key2);
		
		if(iCoeff==-1) {
			a2.put(new JSONObject());
			return a2.getJSONObject(a2.length()-1);
		} else {
			while(a2.length()<iCoeff+1)
				a2.put(new JSONObject());
			return a2.getJSONObject(iCoeff);
		}
	}
	
	private void setAdditionalValues(JSONObject buildStep, JSONObject newElement) {

		Iterator<String> keys = buildStep.keys();

		while( keys.hasNext() ) {
			String key = (String)keys.next();
			if      ( buildStep.get(key) instanceof String     ) newElement.put(key, buildStep.getString(key));
			else if ( buildStep.get(key) instanceof Integer    ) newElement.put(key, buildStep.getInt(key)); 
			else if ( buildStep.get(key) instanceof Boolean    ) newElement.put(key, buildStep.getBoolean(key));
			else if ( buildStep.get(key) instanceof Number     ) newElement.put(key, buildStep.getNumber(key));
			else if ( buildStep.get(key) instanceof JSONObject ) newElement.put(key, ModelComponent.clone_JSON_O(buildStep.getJSONObject(key)));
			else if ( buildStep.get(key) instanceof JSONArray  ) newElement.put(key, ModelComponent.clone_JSON_A(buildStep.getJSONArray(key) ));
			else throw new java.lang.RuntimeException("Can't copy arguments in JSONObject: "+buildStep.toString());
		}
	}

	//Explicitly adds a single count table.
	private void addTable(JSONObject buildStep) {
		
		//System.out.println(">> Running 'addTable'.");
		JSONObject newElement                 = createAndReturnEntry_A("modelSettings",           "countTable", -1);
		createAndReturnEntry_A(                                        "modelFittingConstraints", "countTable", -1);
		
		JSONObject stepCopy                   = ModelComponent.clone_JSON_O(buildStep);
	
		newElement.put( "countTableFile",       stepCopy.getString("countTableFile"));
		if(stepCopy.has("testFolds"))
			newElement.put("testFolds",stepCopy.getJSONArray("testFolds"));
		newElement.put( "variableRegionLength", stepCopy.getInt("variableRegionLength"));
		newElement.put( "nColumns",             stepCopy.getInt("nColumns"));		
		newElement.put( "nColumns",             stepCopy.getInt("nColumns"));
		stepCopy.remove("countTableFile");
		stepCopy.remove("variableRegionLength");
		stepCopy.remove("nColumns");
		
		if(stepCopy.has("function"))
			stepCopy.remove("function");

		setAdditionalValues(stepCopy, newElement);
		
		buildStep.keys();
		
	}

	//Adds a single count table from CELLX
	private void addTableDB(JSONObject buildStep) {

		//Create a new build step using the database and calls addTable to do the work.
		JSONObject newBuildStep = new JSONObject();
		newBuildStep.put("function", "addTable");
		int countTableID     = buildStep.getInt("count_table_id");
		try {
			Connection c     = getDBConnection();

			// GETS INFORMATION ABOUT THE COUNT TABLE
			Statement stmtCT = c.createStatement();
			ResultSet rsCT   = stmtCT.executeQuery("SELECT * FROM count_table WHERE count_table_id="+countTableID+";");
			rsCT.next();
			String  study_name         = rsCT.getString("study_name");
			String  experiment_name    = rsCT.getString("experiment_name");
			Integer max_reads          = rsCT.getInt("max_reads");
			int nColumns               = ((Integer[]) (rsCT.getArray("col_sums")).getArray()).length;
			newBuildStep.put("nColumns", nColumns);
			
			newBuildStep.put("countTableFile", cellx_dir+"/pipeline/countTables/"+study_name+"/"+experiment_name+"."+max_reads+".tsv.gz");
			newBuildStep.put("inputFileType",  "tsv.gz");
			stmtCT.close();
			rsCT.close();

			//GETS INFORMATION ABOUT THE TABLE COLUMNS
			Statement stmtCol = c.createStatement();
			ResultSet rsCol   = stmtCol.executeQuery("SELECT library_name FROM libs_in_expr('"+study_name+"', '"+experiment_name+"');");
			boolean first     = true;
			JSONArray modeledColumns = new JSONArray();
			while( rsCol.next() ) {
				String  library_name = rsCol.getString("library_name");
				
				//GETS INFORMATION ABOUT LIBRARIES. 
				Statement stmtLib = c.createStatement();
				ResultSet rsLib   = stmtLib.executeQuery("SELECT * FROM library WHERE study_name='"+study_name+"' AND library_name='"+library_name+"';");
				rsLib.next();
				modeledColumns.put(rsLib.getInt("round"));
				if(first) {
					newBuildStep.put("variableRegionLength", rsLib.getInt("len"));
					newBuildStep.put("leftFlank",            rsLib.getString("lflank"));
					newBuildStep.put("rightFlank",           rsLib.getString("rflank"));
					
					first = false;
				}
				stmtLib.close();
				rsLib.close();
			}
			rsCol.close();
			
			// GETS INFORMATION ABOUT THE TESTING DATA, IF APPLICABLE
			if(buildStep.has("count_table_id_test")) {
				int countTableIDtest     = buildStep.getInt("count_table_id_test");
				stmtCT = c.createStatement();
				rsCT   = stmtCT.executeQuery("SELECT * FROM count_table WHERE count_table_id="+countTableIDtest+";");
				rsCT.next();
				String  study_name_test      = rsCT.getString("study_name");
				String  experiment_name_test = rsCT.getString("experiment_name");
				Integer max_reads_test       = rsCT.getInt("max_reads");
				newBuildStep.put("countTableFileTest", cellx_dir+"/pipeline/countTables/"+study_name_test+"/"+experiment_name_test+"."+max_reads_test+".tsv.gz");
				stmtCT.close();
				rsCT.close();
			}
				
			//Closes DB connection
			c.close();
			
			if(buildStep.has("modeledColumns"))
				newBuildStep.put("modeledColumns",               buildStep.get("modeledColumns"));
			else
				newBuildStep.put("modeledColumns",               modeledColumns);
			
		} catch ( Exception e ) {
			System.err.println( e.getClass().getName()+": "+ e.getMessage() );
			System.exit(0);
		}
		
		if(buildStep.has("transliterate")) {
			newBuildStep.put("transliterate", buildStep.getJSONObject("transliterate"));
		}

		if(buildStep.has("testFolds")) {
			newBuildStep.put("testFolds", buildStep.getJSONArray("testFolds"));
		}

		addTable(newBuildStep);
	}
	
	//Adds multiple tables forom CellX
	private void addTableDBs(JSONObject buildStep) {
		
		//System.out.println(">> Running 'addTableDBs'.");
		for(int i=0; i<buildStep.getJSONArray("count_table_ids").length(); i++) {
			JSONObject stepCopy = ModelComponent.clone_JSON_O(buildStep);
			stepCopy.put("count_table_id", buildStep.getJSONArray("count_table_ids").getInt(i));
			stepCopy.remove("count_table_ids");
			addTableDB(stepCopy);
		}
	}

	//Explicitly adds a SELEX-experiment enrichment model
	private void addSELEX(JSONObject buildStep) {

		//System.out.println(">> Running 'addSelex'.");
		JSONObject newElement   = createAndReturnEntry_A("modelSettings",           "enrichmentModel", -1);
		createAndReturnEntry_A(                          "modelFittingConstraints", "enrichmentModel", -1);
		JSONObject stepCopy     = ModelComponent.clone_JSON_O(buildStep);
		
		
		//1) Pull modifications from count table!!
		newElement.put("modelType",            buildStep.has("modelType") ? buildStep.get("modelType") : "SELEX");

		//2) Save chemical modifications.
		JSONArray newMod = new JSONArray();
		newElement.put("modifications", newMod);
		if(stepCopy.has("modifications")) {
			JSONArray inMod = stepCopy.getJSONArray("modifications");
			for(int iMod=0; iMod<inMod.length(); iMod++)
				newMod.put(inMod.get(iMod));
		}
		
		//3) Copies binding mode 
		if(stepCopy.has("bindingModes")) {
			newElement.put("bindingModes",            stepCopy.getJSONArray("bindingModes"));
		} else {
			JSONArray newBM  = new JSONArray();
			newBM.put(-1);
			newElement.put("bindingModes",            newBM);
		}
		
		if(stepCopy.has("bindingModeInteractions")) {
			newElement.put("bindingModeInteractions", stepCopy.getJSONArray("bindingModeInteractions"));
		} else {
			JSONArray newInt = new JSONArray();
			newInt.put(-1);
			newElement.put("bindingModeInteractions", newInt);
		}
		
		//Copies round-specific binding modes/interactions
		if(stepCopy.has(    "roundSpecificBindingModes")) {
			newElement.put( "roundSpecificBindingModes",            stepCopy.getJSONArray("roundSpecificBindingModes"));
			stepCopy.remove("roundSpecificBindingModes");
		}
		if(stepCopy.has(    "roundSpecificBindingModeInteractions")) {
			newElement.put( "roundSpecificBindingModeInteractions", stepCopy.getJSONArray("roundSpecificBindingModeInteractions"));
			stepCopy.remove("roundSpecificBindingModeInteractions");
		}

		if(newElement.get("modelType").equals("SELEX_NRLB")) {
			//4) Round is copied for "SELEX_NRLB"
			newElement.put("round", stepCopy.getInt("round"));
			
		} else {
		}

		//5) Write additional arguments, such as "bindingSaturation"
		stepCopy.remove("function");
		stepCopy.remove("bindingModes");
		stepCopy.remove("bindingModeInteractions");
		if(stepCopy.has("modifications"))
			stepCopy.remove("modifications");
		if(stepCopy.has("rounds"))
			stepCopy.remove("rounds");
		setAdditionalValues(stepCopy, newElement);

	}
	
	//Adds a single SELEX experiment using the database.
	private void addSELEXDB(JSONObject buildStep) {

		//Create a new build step using the database and calls addTable to do the work.
		JSONObject newBuildStep = ModelComponent.clone_JSON_O(buildStep);
		newBuildStep.put("function", "addSELEX");
		int countTableID     = buildStep.getInt("count_table_id");

		//We need to add:
		//     "modifications":           { "type": "array",  "items": {"type": "string"}  },
		//     "rounds":                  { "type": "array",  "items": {"type": "integer"} },
		//     "bindingModes":            { "type": "array",  "items": {"type": "integer"}, "default": [-1] },
		//     "bindingModeInteractions": { "type": "array",  "items": {"type": "integer"}, "default": [-1] }
		try {
			Connection c     = getDBConnection();

			// GETS STUDY AND EXPERIMENT NAMES. 
			Statement stmtCT = c.createStatement();
			ResultSet rsCT   = stmtCT.executeQuery("SELECT * FROM count_table WHERE count_table_id="+countTableID+";");
			rsCT.next();
			String  study_name         = rsCT.getString("study_name");
			String  experiment_name    = rsCT.getString("experiment_name");
			stmtCT.close();
			rsCT.close();
			String sqlCommand = "SELECT chem_mod FROM experiment WHERE study_name='"+study_name+"' AND experiment_name='"+experiment_name+"';";

			// CREATES A LIST OF CHEMICAL MODIFICATIONS.
			Statement stmtExp = c.createStatement();
			ResultSet rsExp   = stmtExp.executeQuery(sqlCommand);
			rsExp.next();

			String chem_mod = rsExp.getString("chem_mod");
			JSONArray newMod = new JSONArray();
			newBuildStep.put("modifications", newMod);
			if(chem_mod!=null)
				newMod.put(chem_mod);
			stmtExp.close();
			rsExp.close();

			// CREATES A LIST OF ROUNDS.
			JSONArray newRounds      = new JSONArray();
			Statement stmtCol        = c.createStatement();
			String sqlCommand2       = "SELECT library_name FROM libs_in_expr('"+study_name+"', '"+experiment_name+"');";
			ResultSet rsCol          = stmtCol.executeQuery(sqlCommand2);

			while( rsCol.next() ) {
				String  library_name = rsCol.getString("library_name");
				Statement stmtLib    = c.createStatement();
				ResultSet rsLib      = stmtLib.executeQuery("SELECT * FROM library WHERE study_name='"+study_name+"' AND library_name='"+library_name+"';");
				rsLib.next();
				newRounds.put(rsLib.getInt("round"));
				stmtLib.close();
				rsLib.close();
			}

			if(newBuildStep.has("modelType")&&newBuildStep.getString("modelType").equals("SELEX_NRLB")) {
				newBuildStep.put("round",  newRounds.get(newRounds.length()-1));
			} else {
				newBuildStep.put("rounds", newRounds);
			}

			rsCol.close();
			c.close();

		} catch ( Exception e ) {
			System.err.println( e.getClass().getName()+": "+ e.getMessage() );
			System.exit(0);
		}
		
		if(!newBuildStep.has("bindingModes")) {
			JSONArray newBM  = new JSONArray();
			newBM.put( -1);
			newBuildStep.put("bindingModes",            newBM );
		}
		if(!newBuildStep.has("bindingModeInteractions")) {
			JSONArray newInt = new JSONArray();
			newInt.put(-1);
			newBuildStep.put("bindingModeInteractions", newInt);			
		}

		
		if(newBuildStep.has("count_table_id"))
			newBuildStep.remove("count_table_id");
		addSELEX(newBuildStep);
		
	}

	//Adds both a single count table and a single SELEX experiment.
	private void addSELEXTableDB(JSONObject buildStep) {

		JSONObject stepCopy;
		String[] selexSettings = {"modelType", "modifications", "bindingModes", "bindingModeInteractions", "rounds", "concentration", "bindingSaturation", "r0KUsed", "r0KsTested"};
		String[] countSettings = {"countTableFile", "inputFileType", "nColumns", "variableRegionLength", "rightFlank", "leftFlank", "transliterate"};

		//Makes builds steps for TableDB
		stepCopy = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.put("function", "addTableDB");
		for(int iSet=0; iSet<selexSettings.length; iSet++)
			if(stepCopy.has(selexSettings[iSet]))
				stepCopy.remove(selexSettings[iSet]);

		addTableDB(stepCopy);
		
		//Makes builds steps for SELEXDB
		stepCopy = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.put("function", "addSELEXDB");
		for(int iSet=0; iSet<countSettings.length; iSet++)
			if(stepCopy.has(countSettings[iSet]))
				stepCopy.remove(countSettings[iSet]);
		addSELEXDB(stepCopy);
		

	}
	

	//Adds both multiple count table and a single SELEX experiment.
	private void addSELEXTableDBs(JSONObject buildStep) {		
		
		for(int i=0; i<buildStep.getJSONArray("count_table_ids").length(); i++) {
			JSONObject stepCopy = ModelComponent.clone_JSON_O(buildStep);
			stepCopy.put("count_table_id", buildStep.getJSONArray("count_table_ids").getInt(i));
			stepCopy.remove("count_table_ids");
			
			//Passes the appropriate concentrations to addSELEXTableDB
			if(stepCopy.has("concentrations")) {
				stepCopy.put("concentration", buildStep.getJSONArray("concentrations").getDouble(i));
				stepCopy.remove("concentrations");
			}
			
			//Uses appropriate columns if experiment-specific modeledColumns is specified (of the form [[0,1], [1,3],...])
			if(stepCopy.has("modeledColumns") && stepCopy.getJSONArray("modeledColumns").get(0) instanceof JSONArray )
				stepCopy.put("modeledColumns", stepCopy.getJSONArray("modeledColumns").getJSONArray(i));

			if(stepCopy.has("testFolds") && stepCopy.getJSONArray("testFolds").get(0) instanceof JSONArray )
				stepCopy.put("testFolds", stepCopy.getJSONArray("testFolds").getJSONArray(i));

			addSELEXTableDB(stepCopy);
		}
	}
	
	private void addBindingMode(JSONObject buildStep) {
		
		//System.out.println(">> Running 'addBindingModes'.");
		JSONObject newElement = createAndReturnEntry_A("modelSettings",           "bindingModes", -1);
		createAndReturnEntry_A(                        "modelFittingConstraints", "bindingModes", -1);

		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		
		newElement.put("size",                  buildStep.getInt("size"));

		if(stepCopy.has(    "modifications")) {
			JSONArray newModArray = new JSONArray();
			ModelComponent.clone_JSON_A(newModArray, buildStep.getJSONArray("modifications"));
			newElement.put( "modifications", newModArray);
			stepCopy.remove("modifications");
		}
		
		stepCopy.remove("function");
		stepCopy.remove("size");
		setAdditionalValues(stepCopy, newElement);
	}
	
	//Adds a binding mode with NS binding mode
	private void addNS(JSONObject buildStep) {
		JSONObject newElement  = createAndReturnEntry_A("modelSettings",           "bindingModes", -1);
		createAndReturnEntry_A(                         "modelFittingConstraints", "bindingModes", -1);
		newElement.put("size", 0);
	}
	
	//Adds a binding mode interactions
	private void addInteraction(JSONObject buildStep) {
		
		//System.out.println(">> Running 'addInteraction'.");
		JSONObject newElement  = createAndReturnEntry_A("modelSettings",           "bindingModeInteractions", -1);
		createAndReturnEntry_A(                         "modelFittingConstraints", "bindingModeInteractions", -1);
		JSONObject stepCopy    = ModelComponent.clone_JSON_O(buildStep);
		
		//Saves the binding modes.
		JSONArray BMs = new JSONArray();
		newElement.put("bindingModes", BMs);
		for(int iBM=0; iBM<2; iBM++) 
			BMs.put(stepCopy.getJSONArray("bindingModes").getInt(iBM));
		
		stepCopy.remove("function");
		stepCopy.remove("bindingModes");
		setAdditionalValues(stepCopy, newElement);
		
	}
	
	//Adds output.
	private void output(JSONObject buildStep) {
		JSONObject newElement = createAndReturnEntry_O("optimizerSetting", "output");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		newElement.put( "outputPath",  buildStep.getString("outputPath"));
		stepCopy.remove("outputPath");
		newElement.put( "baseName",    buildStep.getString("baseName"));
		stepCopy.remove("baseName");
		if(stepCopy.has("verbose")) {
			newElement.put("verbose", stepCopy.getBoolean("verbose"));
			stepCopy.remove("verbose");
		} else {
			newElement.put("verbose", true);
		}
		setAdditionalValues(stepCopy, newElement);

	}
	
	//Sets the alphabet
	private void setAlphabet(JSONObject buildStep) {
		if(!config.has("modelSettings"))
			config.put("modelSettings", new JSONObject());
		JSONObject oSett = config.getJSONObject("modelSettings");
		String letterComplement = buildStep.has("letterComplement") ? buildStep.getString("letterComplement") : "C-G,A-T";
		String letterOrder      = buildStep.has("letterOrder")      ? buildStep.getString("letterOrder")      : CombinedLikelihood.alphabetDefToOrderedLetters(letterComplement);
		oSett.put("letterComplement", letterComplement);
		oSett.put("letterOrder",      letterOrder);
	}
	
	
	//Adds output using DB
	private void outputDB(JSONObject buildStep) {
		//System.out.println(">> Running 'outputDB'.");
		JSONObject newBuildStep = ModelComponent.clone_JSON_O(buildStep);
		int fitID     = buildStep.getInt("fit_id");
		newBuildStep.remove("fit_id");
		newBuildStep.put("function", "output");
		
		//We need to add:
		// "outputPath": { "type": "string" }, "baseName": { "type": "string" }
		try {
			Connection c      = getDBConnection();
			Statement stmtFit = c.createStatement();
			ResultSet rsFit   = stmtFit.executeQuery("SELECT * FROM fits WHERE fit_id="+fitID+";");
			rsFit.next();
			newBuildStep.put("outputPath",      cellx_dir+"/pipeline/fits/"+rsFit.getString("fit_path"));
			newBuildStep.put("baseName",        "fit");
			newBuildStep.put("printTrajectory", true);
			newBuildStep.put("verbose", true);
			stmtFit.close();
			rsFit.close();
			c.close();

		} catch ( Exception e ) {
			System.err.println( e.getClass().getName()+": "+ e.getMessage() );
			System.exit(0);
		}


		output(newBuildStep);
		
	}
	
	//Adds optimizer settings.
	private void optimizerSetting(JSONObject buildStep) {

		if(!config.has( "optimizerSetting"))
			config.put("optimizerSetting", new JSONObject());
		JSONObject newElement = config.getJSONObject("optimizerSetting");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		setAdditionalValues(stepCopy, newElement);

	}

	//Adds LBFGS settings.
	private void lbfgsSettings(JSONObject buildStep) {

		JSONObject newElement = createAndReturnEntry_O("optimizerSetting", "lbfgsSettings");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		setAdditionalValues(stepCopy, newElement);
	}
	
	//Adds pattern search settings.
	private void patternSearchSettings(JSONObject buildStep) {

		JSONObject newElement = createAndReturnEntry_O("optimizerSetting", "patternSearchSettings");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		setAdditionalValues(stepCopy, newElement);
	}
	
	
	//Adds S-LBFGS Settings.
	private void slbfgsSettings(JSONObject buildStep) {
		JSONObject newElement = createAndReturnEntry_O("optimizerSetting", "slbfgsSettings");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		setAdditionalValues(stepCopy, newElement);

	}
	
	
	//Adds basic model fitting constraints (i.e. addBindingModesSequentially).
	private void basicConstraints(JSONObject buildStep) {
		
		if(!config.has( "modelFittingConstraints"))
			config.put("modelFittingConstraints", new JSONObject());
		JSONObject newElement = config.getJSONObject("modelFittingConstraints");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		setAdditionalValues(stepCopy, newElement);
	}
	
	//Adds enrichment model constraints.
	private void enrichmentModelConstraints(JSONObject buildStep) {
		
		int index = buildStep.getInt("index");
		if(buildStep.has("index") && buildStep.getInt("index")!=-1) {
			JSONObject newElement = createAndReturnEntry_A("modelFittingConstraints", "enrichmentModel", index);
			buildStep.remove("index");
			JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
			stepCopy.remove("function");
			setAdditionalValues(stepCopy, newElement);
		} else {
			for(int iCurr=0; iCurr<config.getJSONObject("modelSettings").getJSONArray("enrichmentModel").length(); iCurr++) {
				JSONObject oBuildStep = ModelComponent.clone_JSON_O(buildStep);
				oBuildStep.put("index", iCurr);
				enrichmentModelConstraints(oBuildStep);
			}

		}
	}

	//Adds binding mode constraints.
	private void bindingModeConstraints(JSONObject buildStep) {

		JSONObject newElement = createAndReturnEntry_A("modelFittingConstraints", "bindingModes",    buildStep.getInt("index"));
		buildStep.remove("index");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		setAdditionalValues(stepCopy, newElement);
	}

	//Adds interactions constraints.
	private void interactionConstraints(JSONObject buildStep) {

		JSONObject newElement = createAndReturnEntry_A("modelFittingConstraints", "bindingModeInteractions", buildStep.getInt("index"));
		buildStep.remove("index");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		setAdditionalValues(stepCopy, newElement);
	}
	
	private void NRLBConstraints(JSONObject buildStep) {
		//Runs "basicConstraints".
		JSONObject basicSteps = new JSONObject();
		JSONArray flankL      = new JSONArray();
		JSONArray bmSet       = config.getJSONObject("modelSettings").getJSONArray("bindingModes");
		for(int iBM=0; iBM<bmSet.length(); iBM++)
			if(bmSet.getJSONObject(iBM).has("flankLength"))
				flankL.put(bmSet.getJSONObject(iBM).get("flankLength"));
		if(flankL.length()==0)
			flankL.put(0);
		buildStep.put("flankLengths", flankL);
		buildStep.put("nShifts",               buildStep.has("nShifts")               ? buildStep.get("nShifts")               : 0     );
		buildStep.put("singleModeLengthSweep", buildStep.has("singleModeLengthSweep") ? buildStep.get("singleModeLengthSweep") : "true");
		basicConstraints(basicSteps);
			
		//Adds binding mode.
		for(int iBM=0; iBM<bmSet.length(); iBM++) {
			JSONObject bmConst = new JSONObject();
			bmConst.put("function",  "bindingModeConstraints");
			bmConst.put("symmetryString", buildStep.has("symmetryString") ? buildStep.get("symmetryString") : "null");
			if(buildStep.has("startK") || bmSet.getJSONObject(iBM).has("size"))
				bmConst.put("startK", buildStep.has("startK") ? buildStep.get("startK") : Math.max(1,bmSet.getJSONObject(iBM).getInt("size")) );
			if(buildStep.has("endK")   || bmSet.getJSONObject(iBM).has("size"))
				bmConst.put("endK",   buildStep.has("endK")   ? buildStep.get("endK")   : Math.max(1,bmSet.getJSONObject(iBM).getInt("size")) );
			bmConst.put("index", iBM);
			bindingModeConstraints(bmConst);
		}
	}
	
	//Adds binding mode constraints.
	private void symmetry(JSONObject buildStep) {

		int index = buildStep.getInt("index");
		JSONObject newElement = createAndReturnEntry_A("modelFittingConstraints", "bindingModes", index);
		buildStep.remove("index");
		JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
		stepCopy.remove("function");
		if(!stepCopy.has("symmetryString")) {
			int size = config.getJSONObject("modelSettings").getJSONArray("bindingModes").getJSONObject(index).getInt("size");
			stepCopy.put("symmetryString", "1:"+size+":1");
		}
		setAdditionalValues(stepCopy, newElement);
	}
	
	//Adds binding mode seed.
	private void bindingModeSeed(JSONObject buildStep) {

		if(buildStep.has("index")) {
			JSONObject newElement = createAndReturnEntry_A("modelSeeding", "bindingModes",    buildStep.getInt("index"));
			buildStep.remove("index");
			JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
			stepCopy.remove("function");
			setAdditionalValues(stepCopy, newElement);	 
		} else {
			int nBM = config.getJSONObject("modelSettings").getJSONArray("bindingModes").length();
			for(int iBM=0; iBM<nBM; iBM++) {
				JSONObject newElement = createAndReturnEntry_A("modelSeeding", "bindingModes",    iBM);
				JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
				stepCopy.remove("function");
				setAdditionalValues(stepCopy, newElement);	
			}
		}
	}
	
	//Adds seed to enrichment model
	private void enrichmentModelSeed(JSONObject buildStep) {
		
		if(buildStep.has("index") && buildStep.getInt("index")!=-1) {
			JSONObject newElement = createAndReturnEntry_A("modelSeeding", "enrichmentModel",    buildStep.getInt("index"));
			buildStep.remove("index");
			JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
			stepCopy.remove("function");
			setAdditionalValues(stepCopy, newElement);			
		} else {
			int nEnr = config.getJSONObject("modelSettings").getJSONArray("enrichmentModel").length();
			for(int iEnr=0; iEnr<nEnr; iEnr++) {
				JSONObject newElement = createAndReturnEntry_A("modelSeeding", "enrichmentModel",  iEnr);
				JSONObject stepCopy   = ModelComponent.clone_JSON_O(buildStep);
				stepCopy.remove("function");
				setAdditionalValues(stepCopy, newElement);			
			}
		}
	}

}

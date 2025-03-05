from features import calculate_features
import numpy as np
import pandas as pd
import joblib
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def predict_compound(smiles, test):
    try:
        model = joblib.load(f"models/{test}.lb")
        features = calculate_features(smiles)
        
        if features is None:
            return {"error": "Invalid SMILES string."}
        
        # Reindex to match model columns
        model_columns = model.feature_names_in_
        features = features.reindex(columns=model_columns, fill_value=0)
        
        # Convert to numpy array and ensure 2D shape (1, n_features)
        feature_array = features.to_numpy()
        if feature_array.ndim > 2:
            feature_array = feature_array.reshape(1, -1)
        elif feature_array.ndim == 1:
            feature_array = feature_array.reshape(1, -1)
        
        # Ensure final shape is (1, n_features)
        if feature_array.shape[0] != 1 or feature_array.ndim != 2:
            raise ValueError(f"Feature array shape {feature_array.shape} is not 2D with one row.")
        
        prediction = model.predict(feature_array)
        
        # Convert all features to a dictionary with native Python types
        result = {col: features[col].iloc[0].item() for col in features.columns}
        result["SMILES"] = smiles
        result["prediction"] = "ACTIVE" if prediction[0] == 1 else "INACTIVE"
        result["test"] = test
        
        # Calculate drug-likeness (Lipinski's Rule of Five)
        violations = 0
        if result['XLogP'] > 5: violations += 1
        if result['MW'] > 500: violations += 1
        if result['NumHBA'] > 10: violations += 1
        if result['NumHBD'] > 5: violations += 1
        result["Lipinski_Violations"] = violations
        result["Druglike"] = "Yes" if violations <= 1 else "No"
        
        # Generate 2D molecule image with SMILES as label, increased height for label spacing
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol, size=(400, 300), legend=smiles)  # Increased height to 400
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode('utf-8')
        result["molecule_image"] = f"data:image/png;base64,{img_str}"
        
        return result
        
    except FileNotFoundError:
        return {"error": f"Model file 'models/{test}.lb' not found."}
    except Exception as e:
        return {"error": str(e)}

if __name__ == "__main__":
    example_smiles = "CN1CCC2=C(C1)[C@H](C3=CC(=CC=C3F)OC2)OCO"
    print(predict_compound(example_smiles, 'AID456'))
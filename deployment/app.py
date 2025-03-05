from flask import Flask, request, jsonify, render_template
import os

app = Flask(__name__, static_folder='static', template_folder='templates')

# Mock ML prediction function 
def predict_compound(smiles, test):
    is_active = len(smiles) > 10
    
    # Mock molecular properties 
    return {
        "XLogP": 3.45,
        "PSA": 78.9,
        "NumRot": 5,
        "NumHBA": 3,
        "NumHBD": 2,
        "MW": 324.5,
        "prediction": "ACTIVE" if is_active else "INACTIVE"
    }

# Route to serve the frontend
@app.route('/')
def serve_frontend():
    return render_template('index.html')

# API endpoint for prediction
@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    smiles = data.get('smiles', '')
    test = data.get('test', '')

    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    # Get prediction and properties
    result = predict_compound(smiles, test)
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True)
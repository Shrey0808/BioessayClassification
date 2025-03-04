from flask import Flask, render_template, request, jsonify
from prediction import predict_compound

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    smiles = data.get('smiles')
    test = data.get('test')
    
    if not smiles or not test:
        return jsonify({"error": "SMILES and test type are required."}), 400
    
    result = predict_compound(smiles, test)
    
    if "error" in result:
        return jsonify(result), 400
    
    return jsonify(result), 200

if __name__ == '__main__':
    app.run(debug=True)
from flask import Flask, request, jsonify, render_template, send_file
from prediction import predict_compound
import csv
from io import StringIO, BytesIO

app = Flask(__name__)

# Store the latest result globally
latest_result = None

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    global latest_result
    data = request.get_json()
    smiles = data['smiles']
    test = data['test']
    result = predict_compound(smiles, test)
    latest_result = result if 'error' not in result else None  # Store successful results
    # Return subset for UI display
    if latest_result:
        ui_result = {
            "XLogP": latest_result["XLogP"],
            "PSA": latest_result["PSA"],
            "NumRot": latest_result["NumRot"],
            "NumHBA": latest_result["NumHBA"],
            "NumHBD": latest_result["NumHBD"],
            "MW": latest_result["MW"],
            "prediction": latest_result["prediction"],
            "Lipinski_Violations": latest_result["Lipinski_Violations"],
            "Druglike": latest_result["Druglike"],
            "molecule_image": latest_result["molecule_image"]
        }
        return jsonify(ui_result)
    return jsonify(result)

@app.route('/download', methods=['GET'])
def download():
    global latest_result
    if latest_result is None:
        return jsonify({"error": "No analysis result available to download."}), 400
    
    # Create a CSV file with all features
    output = StringIO()
    # Exclude molecule_image from CSV as it's not tabular data
    csv_result = {k: v for k, v in latest_result.items() if k != "molecule_image"}
    writer = csv.DictWriter(output, fieldnames=csv_result.keys())
    writer.writeheader()
    writer.writerow(csv_result)
    
    buffer = BytesIO()
    buffer.write(output.getvalue().encode('utf-8'))
    buffer.seek(0)
    output.close()
    
    return send_file(
        buffer,
        as_attachment=True,
        download_name='analysis_result.csv',
        mimetype='text/csv'
    )

if __name__ == '__main__':
    app.run(debug=True)
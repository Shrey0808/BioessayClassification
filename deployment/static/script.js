window.selectedValue = 'AID362';
let latestPredictionData = null;

async function generatePrediction() {
    console.log('Analyze button clicked');
    const smiles = document.querySelector('.smiles-input').value;
    const test = window.selectedValue;
    const loader = document.getElementById('loader');
    const results = document.getElementById('results');
    const generateBtn = document.querySelector('.generate-btn');

    if (!smiles) {
        console.log('No SMILES input provided');
        document.body.classList.add('popup-active');
        await Swal.fire({
            icon: 'warning',
            title: 'Input Required',
            text: 'Please provide a valid SMILES notation to proceed.',
            background: '#D4E8F8',
            color: '#3282B8',
            confirmButtonColor: '#3282B8',
            customClass: {
                popup: 'rounded-popup'
            }
        });
        document.body.classList.remove('popup-active');
        return;
    }

    generateBtn.classList.add('loading');
    loader.style.display = 'block';
    results.style.display = 'none';

    try {
        console.log('Sending request to /predict:', { smiles, test });
        const response = await fetch('/predict', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles, test })
        });

        const data = await response.json();
        console.log('Response received:', data);

        generateBtn.classList.remove('loading');
        loader.style.display = 'none';

        if (!response.ok || data.error) {
            console.log('Error detected:', data.error);
            document.body.classList.add('popup-active');
            const errorMessage = data.error === 'Invalid SMILES string.' 
                ? 'The input SMILES string is invalid.' 
                : data.error || 'An unexpected error occurred.';
            await Swal.fire({
                icon: 'error',
                title: 'Error',
                text: errorMessage,
                background: '#D4E8F8',
                color: '#3282B8',
                confirmButtonColor: '#3282B8',
                customClass: {
                    popup: 'rounded-popup'
                }
            });
            document.body.classList.remove('popup-active');
            return;
        }

        console.log('Success case triggered');
        latestPredictionData = data; // Store for download
        document.body.classList.add('popup-active');
        const result = await Swal.fire({
            icon: 'success',
            title: 'Analysis Complete',
            text: 'The compound has been analyzed successfully.',
            background: '#D4E8F8',
            color: '#3282B8',
            confirmButtonColor: '#3282B8',
            customClass: {
                popup: 'rounded-popup'
            }
        });

        if (result.isConfirmed) {
            results.style.display = 'block';
            
            // Display molecule image
            const moleculeImage = document.getElementById('molecule-image');
            moleculeImage.innerHTML = `<img src="${data.molecule_image}" alt="Molecule Structure" style="max-width: 200px;">`;

            // Compound analysis
            const compoundData = document.getElementById('compound-data');
            compoundData.innerHTML = `
                <span>XLogP: ${data.XLogP}</span>
                <span>PSA: ${data.PSA} Å²</span>
                <span>Rotatable Bonds: ${data.NumRot}</span>
                <span>H-Bond Acceptors: ${data.NumHBA}</span>
                <span>H-Bond Donors: ${data.NumHBD}</span>
                <span>Molecular Weight: ${data.MW} g/mol</span>
            `;
            const spans = compoundData.querySelectorAll('span');
            spans.forEach((span, index) => {
                span.style.animation = `fadeSlide 0.5s ease-out ${index * 0.2}s forwards`;
            });

            // Drug-likeness
            const druglikenessData = document.getElementById('druglikeness-data');
            druglikenessData.innerHTML = `
                <span>Lipinski Violations: ${data.Lipinski_Violations}</span>
                <span>Druglike: ${data.Druglike}</span>
            `;
            const drugSpans = druglikenessData.querySelectorAll('span');
            drugSpans.forEach((span, index) => {
                span.style.animation = `fadeSlide 0.5s ease-out ${index * 0.2}s forwards`;
            });

            // Prediction
            const predictionResult = document.getElementById('prediction-result');
            predictionResult.textContent = data.prediction === 'ACTIVE' ? 'Active' : 'Inactive';
            predictionResult.className = 'prediction-result ' + (data.prediction === 'ACTIVE' ? 'active' : 'inactive');
            predictionResult.style.animation = 'cardFlip 0.8s ease-out, glitch 0.3s ease 0.8s';

            document.getElementById('download-btn').style.display = 'block';
        }
        document.body.classList.remove('popup-active');

    } catch (error) {
        console.error('Fetch error:', error);
        generateBtn.classList.remove('loading');
        loader.style.display = 'none';
        document.body.classList.add('popup-active');
        await Swal.fire({
            icon: 'error',
            title: 'Error',
            text: error.message || 'An unexpected error occurred.',
            background: '#D4E8F8',
            color: '#3282B8',
            confirmButtonColor: '#3282B8',
            customClass: {
                popup: 'rounded-popup'
            }
        });
        document.body.classList.remove('popup-active');
    }
}

document.getElementById('download-btn').addEventListener('click', async () => {
    console.log('Download button clicked');
    if (!latestPredictionData) {
        document.body.classList.add('popup-active');
        await Swal.fire({
            icon: 'warning',
            title: 'No Data',
            text: 'No analysis result available to download.',
            background: '#D4E8F8',
            color: '#3282B8',
            confirmButtonColor: '#3282B8',
            customClass: {
                popup: 'rounded-popup'
            }
        });
        document.body.classList.remove('popup-active');
        return;
    }

    try {
        const response = await fetch('/download', {
            method: 'GET'
        });

        if (!response.ok) {
            throw new Error('Failed to download the file.');
        }

        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = 'analysis_result.csv';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        window.URL.revokeObjectURL(url);

        console.log('Download initiated successfully');
    } catch (error) {
        console.error('Download error:', error);
        document.body.classList.add('popup-active');
        await Swal.fire({
            icon: 'error',
            title: 'Download Failed',
            text: error.message || 'An error occurred while downloading the file.',
            background: '#D4E8F8',
            color: '#3282B8',
            confirmButtonColor: '#3282B8',
            customClass: {
                popup: 'rounded-popup'
            }
        });
        document.body.classList.remove('popup-active');
    }
});
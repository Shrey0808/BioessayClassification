window.selectedValue = 'test1'; // Global variable for dropdown selection

async function generatePrediction() {
    const smiles = document.querySelector('.smiles-input').value;
    const test = window.selectedValue; // Use global selected value
    const loader = document.getElementById('loader');
    const results = document.getElementById('results');
    const generateBtn = document.querySelector('.generate-btn');

    if (!smiles) {
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
        const response = await fetch('/predict', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles, test })
        });

        const data = await response.json();

        generateBtn.classList.remove('loading');
        loader.style.display = 'none';

        if (response.ok) {
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

                const predictionResult = document.getElementById('prediction-result');
                predictionResult.textContent = data.prediction === 'ACTIVE' ? 'Active' : 'Inactive';
                predictionResult.className = 'prediction-result ' + (data.prediction === 'ACTIVE' ? 'active' : 'inactive');
                predictionResult.style.animation = 'cardFlip 0.8s ease-out, glitch 0.3s ease 0.8s';

                document.getElementById('download-btn').style.display = 'block';
            }
            document.body.classList.remove('popup-active');
        } else {
            throw new Error(data.error || 'An unexpected error occurred.');
        }
    } catch (error) {
        generateBtn.classList.remove('loading');
        loader.style.display = 'none';

        document.body.classList.add('popup-active');
        await Swal.fire({
            icon: 'error',
            title: 'Error',
            text: error.message,
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
    document.body.classList.add('popup-active');
    await Swal.fire({
        icon: 'info',
        title: 'Download Unavailable',
        text: 'This feature is under development.',
        background: '#D4E8F8',
        color: '#3282B8',
        confirmButtonColor: '#3282B8',
        customClass: {
            popup: 'rounded-popup'
        }
    });
    document.body.classList.remove('popup-active');
});
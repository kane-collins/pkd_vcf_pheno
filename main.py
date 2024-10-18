import os
import logging
import pandas as pd
import json
import io
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from flask import Flask, render_template, request, redirect, url_for, Response, send_file, flash
import plotly
import plotly.express as px
from werkzeug.utils import secure_filename
from werkzeug.security import generate_password_hash, check_password_hash

from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import DataRequired
from flask_login import LoginManager, UserMixin


class LoginForm(FlaskForm):
    username = StringField('Username', validators=[DataRequired()])
    password = StringField('Password', validators=[DataRequired()])
    remember_me = BooleanField('Remember Me')
    submit = SubmitField('Sign In')


# Set matplotlib to use 'Agg' backend for environments without a display
matplotlib.use('Agg')

class Config:
    UPLOAD_FOLDER = os.getenv('UPLOAD_FOLDER', 'uploads')
    PEDIGREE_FOLDER = os.getenv('PEDIGREE_FOLDER', os.path.join(UPLOAD_FOLDER, 'pedigrees'))
    VCF_FOLDER = os.getenv('VCF_FOLDER', os.path.join(UPLOAD_FOLDER, 'vcfs'))
    PARTICIPANT_DATA_FILE = os.getenv('PARTICIPANT_DATA_FILE', 'beacon_format/pkd_individuals.csv')
    #PARTICIPANT_DATA_FILE = os.getenv('PARTICIPANT_DATA_FILE', 'beacon_format/individuals.csv')

    #For user login
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'

# class User(UserMixin): # also db.Model
#     # ...

#     def set_password(self, password):
#         self.password_hash = generate_password_hash(password)

#     def check_password(self, password):
#         return check_password_hash(self.password_hash, password)


app = Flask(__name__)
app.config.from_object(Config)
#login = LoginManager(app)

# Initialize logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def load_participant_data(filepath):
    try:
        df = pd.read_csv(filepath)
        logging.info(f"Loaded participant data from {filepath}")
        return df
    except FileNotFoundError:
        logging.error(f"File not found: {filepath}")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        logging.error(f"No data in file: {filepath}")
        return pd.DataFrame()
    except Exception as e:
        logging.error(f"Error loading participant data: {e}")
        return pd.DataFrame()

# Load the participant data CSV into a DataFrame
input_file_path = app.config['PARTICIPANT_DATA_FILE']
df = load_participant_data(input_file_path)



@app.route('/', methods=['GET', 'POST'])
def index():
    try:
        vcf_folder = app.config['VCF_FOLDER']
        vcf_count = len([f for f in os.listdir(vcf_folder) if f.endswith('.vcf')])
        unique_id_count = df['id'].nunique()
    except Exception as e:
        logging.error(f"Error in index route: {e}")
        vcf_count = unique_id_count = 0

    return render_template('index.html', vcf_count=vcf_count, unique_id_count=unique_id_count)



#Note that this login still isn't working, need to work on this
@app.route('/login', methods=['GET', 'POST'])
def login():
   form = LoginForm()
   if form.validate_on_submit():
       flash('Login requested for user {}, remember_me={}'.format(form.username.data, form.remember_me.data))
       return redirect(url_for('index'))
   return render_template('login.html', title = "Sign In", form = form)




@app.route('/view_phenotype', methods=["GET", "POST"])
def view_phenotype():

    return render_template('view_phenotype.html', column_names=df.columns.values.tolist(), row_data=df.values.tolist())

# Route to search by individual ID
@app.route('/search/<search_id>')
def search(search_id):
    filtered_df = df[df['id'] == search_id].fillna('').loc[:, (df != '').any(axis=0)]
    result = filtered_df.to_dict(orient='records')

    #check if the individual exists
    if result is None:
        return "Individual not found.", 404

    #return the parent IDs, I also want to get the IDs of all family members, along with their relationships, this will be a bit more tricky, will need some sort of recursive search
    mother_id = filtered_df['mother_ID']
    father_id = filtered_df['father_ID']

    family_id = filtered_df['FID']


    vcf_file_path = app.config['VCF_FOLDER']

    filepath=vcf_file_path + '/' + search_id + '.vcf'
    print(filepath)

    # Check if the VCF file exists
    if not os.path.exists(filepath):
        return f"No VCF file available for individual {search_id}.", 404

    var_table = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt'])
    with open(filepath, 'r') as file:
        for line in file: 
            if line.startswith('##') or line.startswith('#'):
                continue
            columns = line.strip().split('\t')
            new_row = pd.DataFrame({'chrom': [columns[0]], 'pos': [columns[1]], 'ref': [columns[3]], 'alt': [columns[4]]})
            var_table = pd.concat([var_table, new_row])
    new_var_table = var_table.to_dict(orient='records')    
    #Need to modify this to work if there is no vcf file

    
    return new_var_table, result, search_id, mother_id, father_id, family_id

@app.route('/clickable_participant_search/<search_id>')
def clickable_participant_search(search_id):
    new_var_table, result, search_id, mother_id, father_id, family_id = search(search_id)
    return render_template('participant_search_v2.html', new_var_table=new_var_table, result=result, search_id=search_id, mother_id=mother_id, father_id=father_id, family_id=family_id)



@app.route('/participant_search', methods=['GET', 'POST'])
def participant_search():
    if request.method == 'POST':
        search_id = request.form['search_id']
        new_var_table, result, search_id, mother_id, father_id, family_id = search(search_id)
        return render_template('participant_search_v2.html', new_var_table=new_var_table, result=result, search_id=search_id, mother_id=mother_id, father_id=father_id, family_id=family_id)

    return render_template('participant_search.html')



@app.route('/variant_search', methods=['GET', 'POST'])
def variant_search():
    result = None
    if request.method == 'POST':
        chrom = request.form['chrom']
        pos = request.form['pos']
        ref = request.form['ref']
        alt = request.form['alt']
        vcf_folder = app.config['VCF_FOLDER']
        result = search_variant(vcf_folder, chrom, pos, ref, alt)
    return render_template('variant_search.html', result=result)




def search_variant(vcf_folder, chrom, pos, ref, alt):
    variant_count = 0
    files_with_variant = 0
    total_files = 0
    variant_files = []

    for filename in os.listdir(vcf_folder):
        if filename.endswith(".vcf"):
            total_files += 1
            filepath = os.path.join(vcf_folder, filename)
            with open(filepath, 'r') as file:
                found_in_file = False
                for line in file:
                    if line.startswith("##") or line.startswith("#"):
                        continue
                    columns = line.strip().split('\t')
                    if columns[0] == chrom and columns[1] == pos and columns[3] == ref and columns[4] == alt:
                        variant_count += 1
                        found_in_file = True
                if found_in_file:
                    files_with_variant += 1
                    variant_files.append(filename)

    percentage = (files_with_variant / total_files) * 100 if total_files > 0 else 0

    gnomad_web_address = 'https://gnomad.broadinstitute.org/variant/' + chrom + '-' + pos + '-' + ref + '-' + alt +'?dataset=gnomad_r4'

    return variant_count, files_with_variant, total_files, percentage, variant_files, gnomad_web_address


ALLOWED_EXTENSIONS = {'png', 'jpg', 'jpeg', 'gif'}
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # Max 16MB upload size

# Helper function to check file extensions
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/add_participant', methods=['GET', 'POST'])
def add_participant():
    global df  # Declare df as global at the beginning of the function

    if request.method == 'POST':
        # Get the file upload
        if 'image' not in request.files:
            return 'No file part'
        file = request.files['image']

        # If the user does not select a file
        if file.filename == '':
            return 'No selected file'

        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], 'imaging', filename))



        # Extract data from the form
        new_data = {}
        for field in request.form:
            new_data[field] = request.form[field]

        # Auto-increment ID
        if not df.empty:
            # Ensure 'id' column is treated as strings and fill NaN with empty string
            df['id'] = df['id'].fillna('').astype(str)
            # Use regex to find the numeric part of the ID
            df['numeric_part'] = df['id'].str.extract('(\d+)', expand=False).fillna('0').astype(int)
            last_id_numeric = df['numeric_part'].max()
            new_id = f"GOI_{last_id_numeric + 1}"
        else:
            new_id = "GOI_1"

        new_data['id'] = new_id

        # Create a DataFrame from the form data
        new_df = pd.DataFrame([new_data])

        # Append the new data to the existing DataFrame
        df = pd.concat([df, new_df], ignore_index=True)

        # Save the updated DataFrame to the CSV file
        try:
            df.to_csv(app.config['PARTICIPANT_DATA_FILE'], index=False)
            logging.info('Data successfully appended to the CSV file')
        except Exception as e:
            logging.error(f"Error saving data to CSV file: {e}")

        return redirect(url_for('index'))
    return render_template('add_participant.html')




# # Configure the folder where your images are stored (outside of static)
# app.config['PEDIGREE_FOLDER'] = '/uploads/pedigrees'

# @app.route('/display/<family_id>')
# def display_image(family_id):
#     # Construct the image URL based on a custom route
#     image_url = f'/pedigree_image/{family_id}.png'
#     return render_template('participant_search_v2.html', image_url=image_url)

# # Route to serve the image directly from the custom folder
# @app.route('/pedigree_image/<filename>')
# def pedigree_image(filename):
#     # Serve the image file from the custom folder
#     return send_from_directory(app.config['PEDIGREE_FOLDER'], filename)


# @app.route('/uploads/pedigrees/<family_id>.png')
# def serve_image(family_id):
#     try:
#         ped_file_path = app.config['PEDIGREE_FOLDER']
#         # Use send_from_directory to serve the image from the uploads folder
#         return send_from_directory(ped_file_path, f"{family_id}.png")
#     except FileNotFoundError:
#         # If the image file doesn't exist, return a 404 error
#         abort(404)

# Route to render the template and pass the family_id to it
@app.route('/display/<family_id>')
def display_image(family_id):
    ped_file_path = app.config['PEDIGREE_FOLDER']

    # Here we construct the image URL, but now it refers to the new route
    image_url = ped_file_path + f'uploads/pedigrees/{family_id}.png'
    print(image_url)
    return render_template('participant_search_v2.html', image_url=image_url)



# @app.route('/display/<family_id>')
# def display_image(family_id):
#     # The image file path is dynamically created
#     image_path = f'uploads/pedigrees/{family_id}.png'
#     return render_template('participant_search_v2.html', image_path=image_path)



# #For the moment, let's just try to add manual pedigree files
# @app.route('/plot_pedigree.png')
# def pedigree_plot(family_id):
#     ped_file_path = app.config['PEDIGREE_FOLDER']

#     filepath=ped_file_path + '/' + family_id + '.vcf'
#     print(filepath)

#     # Check if the VCF file exists
#     if not os.path.exists(filepath):
#         return f"No VCF file available for individual {search_id}.", 404
#     pedigree_plot = open(filepath, 'r') 
#     return render_template("particpant_search_v2.html", pedigree_plot)



@app.route('/plot.png')
def generate_plot():

    fig, axs = plt.subplots(2, 2)  # Create a 2x2 grid of subplots

    # Plot something in each subplot
    counts = df['Sex'].value_counts()
    axs[0, 0].bar(counts.index, counts.values)
    axs[0, 0].set_title("Frequency of male vs female")


    counts = df['KF'].value_counts()
    axs[0, 1].bar(counts.index, counts.values)
    axs[0, 1].set_title("Frequency of kidney failure")


    counts = df['Diag_var'].value_counts()
    axs[1, 0].bar(counts.index, counts.values)
    axs[1, 0].set_title("Diagnostic variant types")


    axs[1, 1].hist(df['age_KF'])
    axs[1, 1].set_title("Age at kidney failure")


    # Tight layout for better spacing
    plt.tight_layout()

    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)  # Go to the start of the file

    return send_file(img, mimetype='image/png')



if __name__ == '__main__':
    if not os.path.exists(app.config['UPLOAD_FOLDER']):
        os.makedirs(app.config['UPLOAD_FOLDER'])
    if not os.path.exists(app.config['VCF_FOLDER']):
        os.makedirs(app.config['VCF_FOLDER'])
    app.run(debug=True, host='127.0.0.1', port=5000)

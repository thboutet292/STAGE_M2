#!/usr/bin/env python3
import os
import re
import shutil
import subprocess
import tempfile

from flask import (
    Flask, request, render_template,
    redirect, url_for, Response,
    stream_with_context, send_from_directory, abort
)

BASE_DIR   = os.path.dirname(os.path.dirname(__file__))
BIN_DIR    = os.path.join(BASE_DIR, 'bin')
DATA_DIR   = os.path.join(BASE_DIR, 'data')
RESULT_DIR = os.path.join(BASE_DIR, 'result')
os.makedirs(RESULT_DIR, exist_ok=True)

app = Flask(__name__, template_folder=os.path.join(BASE_DIR, 'templates'))

def safe(fn):
    return os.path.basename(fn)

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/run', methods=['POST'])
def run():
    # 1) store uploads
    work = tempfile.mkdtemp(prefix='scds_')
    fasta = request.files['fasta']
    fasta_path = os.path.join(work, safe(fasta.filename))
    fasta.save(fasta_path)

    embl_dir = os.path.join(work, 'embl')
    os.makedirs(embl_dir, exist_ok=True)
    for f in request.files.getlist('embl_files'):
        f.save(os.path.join(embl_dir, safe(f.filename)))

    # 2) user options
    outname = safe(request.form['output_name'])
    use_up  = request.form.get('use_updated') == 'on'

    # 3) build pipeline command
    cmd = [
        'bash',
        os.path.join(BIN_DIR, 'pipeline_microannot_sCDS.sh'),
        fasta_path,
        embl_dir,
        outname
    ]
    if use_up:
        cmd.append('--use-updated')

    # number of “Step N:” lines your script prints
    TOTAL_STEPS = 7

    def generate():
        # kickoff
        yield "<script>parent.updateProgress(0,'Initializing…');</script>\n"
        proc = subprocess.Popen(
            cmd, cwd=BASE_DIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True, bufsize=1
        )

        # stream stdout→browser
        for line in proc.stdout:
            esc = line.rstrip().replace("\\", "\\\\").replace("'", "\\'")
            yield f"<script>parent.appendLog('{esc}\\n');</script>\n"

            m = re.match(r'Step\s+(\d+):\s*(.+)', line)
            if m:
                idx, msg = int(m.group(1)), m.group(2).strip()
                pct = int(idx / TOTAL_STEPS * 100)
                yield f"<script>parent.updateProgress({pct},'{msg}');</script>\n"

        proc.wait()

        # copy the final output into result/
        src = os.path.join(BASE_DIR, outname)
        if proc.returncode == 0 and os.path.exists(src):
            dst = os.path.join(RESULT_DIR, outname)
            shutil.copy(src, dst)
            link = url_for('download', filename=outname)
            yield f"<script>parent.finishDownload('{link}');</script>\n"
        else:
            yield "<script>parent.appendLog('⚠️ no output generated');</script>\n"

    return Response(stream_with_context(generate()), mimetype='text/html')


@app.route('/download/<path:filename>')
def download(filename):
    # serve from result/
    return send_from_directory(RESULT_DIR, filename, as_attachment=True)


@app.route('/databases')
def databases():
    originals = []
    # original ortholog DB
    orth_dir = os.path.join(DATA_DIR, 'database_ortholog')
    for f in sorted(os.listdir(orth_dir)):
        originals.append((
            'Ortholog DB', f,
            url_for('download_db', category='ortholog', filename=f)
        ))
    # original non-ortholog DB
    non_dir = os.path.join(DATA_DIR, 'database_not_ortholog')
    for f in sorted(os.listdir(non_dir)):
        originals.append((
            'Non-ortholog DB', f,
            url_for('download_db', category='nonortholog', filename=f)
        ))

    # any generated DBs in result/
    generated = []
    for f in sorted(os.listdir(RESULT_DIR)):
        generated.append((
            f,
            url_for('download_db', category='generated', filename=f)
        ))

    return render_template('databases.html',
                           originals=originals,
                           generated=generated)


@app.route('/download_db/<category>/<path:filename>')
def download_db(category, filename):
    if category == 'ortholog':
        folder = os.path.join(DATA_DIR, 'database_ortholog')
    elif category == 'nonortholog':
        folder = os.path.join(DATA_DIR, 'database_not_ortholog')
    elif category == 'generated':
        folder = RESULT_DIR
    else:
        abort(404)
    return send_from_directory(folder, filename, as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True, port=5001)

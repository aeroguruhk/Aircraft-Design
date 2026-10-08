"""Browser adapter. Runs the same desktop model in an isolated virtual directory.
No uploaded code is executed; the input is JSON only. No files leave the browser.
"""
import base64
import contextlib
import io
import json
import math
from pathlib import Path
import shutil
import sys
import traceback
import zipfile

MODEL=Path('model').resolve()
sys.path.insert(0,str(MODEL))
from run_sizing import run


def validate_browser_config(config):
    if not isinstance(config,dict):
        raise ValueError('The configuration must be a JSON object')
    def finite(value):
        if isinstance(value,(int,float)) and (not math.isfinite(value) or abs(value)>1e9):
            raise ValueError('All numerical inputs must be finite and within teaching-model limits')
        if isinstance(value,dict):
            for item in value.values(): finite(item)
        elif isinstance(value,list):
            for item in value: finite(item)
    finite(config)
    for key in ('aircraft','mission','ws_range'):
        if not isinstance(config.get(key),dict):
            raise ValueError(f'Missing or invalid {key} section')
    if not isinstance(config.get('cases'),list) or not 1<=len(config['cases'])<=20:
        raise ValueError('Provide between 1 and 20 constraints')
    grid=config['ws_range']
    step=grid.get('step_pa',0);start=grid.get('start_pa',0);end=grid.get('end_pa',0)
    if step<=0 or not 0<start<end or (end-start)/step>500:
        raise ValueError('Wing loading range requires a positive step and at most 501 stations')
    a,m=config['aircraft'],config['mission']
    if not 0<=m.get('Range',0)<=10000 or not 0<=m.get('Endurance',0)<=24:
        raise ValueError('Browser examples support 0–10,000 km range and 0–24 hours loiter')
    if not 0<a.get('V_h',0)<=5 or not 0<a.get('V_v',0)<=1:
        raise ValueError('Tail volume coefficients must be positive (V_h ≤ 5; V_v ≤ 1)')
    if not 0<=a.get('cruise_altitude_m',0)<=15000:
        raise ValueError('Cruise altitude must be within 0–15,000 m')


def run_model(config_json):
    log=io.StringIO()
    try:
        if len(config_json)>150000:
            raise ValueError('JSON configuration is too large')
        config=json.loads(config_json)
        validate_browser_config(config)
        work=Path('browser_work').resolve()
        if work.exists(): shutil.rmtree(work)
        work.mkdir()
        path=work/'input.json';path.write_text(json.dumps(config),encoding='utf-8')
        with contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
            summary=run(path,work/'results')
        plots={}
        for key,filename in [('constraints','constraint_diagram.png'),('wing','wing_loads.png'),
                             ('performance','performance_altitude.png'),('range','cruise_range.png')]:
            plots[key]='data:image/png;base64,'+base64.b64encode((work/'results'/filename).read_bytes()).decode('ascii')
        stream=io.BytesIO()
        with zipfile.ZipFile(stream,'w',zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('input.json',json.dumps(config,indent=2))
            for item in sorted((work/'results').iterdir()):
                if item.is_file(): archive.write(item,'results/'+item.name)
            archive.writestr('ABOUT.txt','Aircraft Design Lab — browser run\nAll calculations use the included teaching model.\nSee the online model guide for equations, units and assumptions.\n')
        return json.dumps({'ok':True,'summary':summary,'plots':plots,
                           'zip_base64':base64.b64encode(stream.getvalue()).decode('ascii'),
                           'log':log.getvalue()},allow_nan=False)
    except Exception as exc:
        return json.dumps({'ok':False,'error':str(exc),'log':log.getvalue()+'\n'+traceback.format_exc()})

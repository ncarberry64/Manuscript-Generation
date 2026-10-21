"""Saved early optical means on matching native65x65 history domains; no refit."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'code'))
import run_he0435_joint_scene as prior
js=prior.js
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021/common_early/pilot_means'


def main():
    OUT.mkdir(exist_ok=True);records,fields,eta0,wcss,origin=prior.load_inputs()
    exposures=sorted(wcss);recordmap={(r['exposure'],r['patch']):r for r in records};fieldmap={r['exposure']:r for r in fields}
    anchors=[]
    for r in records[:4]:
        anchors.append(prior.tangent_plane_arcsec(np.asarray(wcss[r['exposure']](*r['guess'])[:2]),origin))
    reports={}
    for mode in ['all_data','training_plus_held']:
        if mode=='all_data':
            solution=np.load(OUT.parent/'all_data/solution.npz');target=js.SceneProblem(records);p=solution['target_parameters'];jointp=solution['parameters'];eta=solution['eta'];fieldstart=target.size+len(eta)
            targets={r['name']:p[3*i:3*i+3] for i,r in enumerate(records)}
            sky=target.sky(p);backgrounds={e:p[target.bg_start+3*i:target.bg_start+3*i+3] for i,e in enumerate(target.exposures)}
        else:
            trainrecords=[r for r in records if r['exposure']!=exposures[-1]];solution=np.load(OUT.parent/'training/solution.npz');target=js.SceneProblem(trainrecords);p=solution['target_parameters'];jointp=solution['parameters'];eta=solution['eta'];fieldstart=target.size+len(eta)
            targets={r['name']:p[3*i:3*i+3] for i,r in enumerate(trainrecords)};sky=target.sky(p);backgrounds={e:p[target.bg_start+3*i:target.bg_start+3*i+3] for i,e in enumerate(target.exposures)}
            held=np.load(OUT.parent/'training/held_solution.npz');hp=held['parameters'];hr=[r for r in records if r['exposure']==exposures[-1]]
            targets.update({r['name']:hp[3*i:3*i+3] for i,r in enumerate(hr)});backgrounds[exposures[-1]]=hp[12:15]
        mode_reports=[]
        for e in exposures:
            historypath=ROOT/'data/strong_lensing_response_completion_20261021/histories'/e.replace('_cal.fits','_histories.npz');history=np.load(historypath)
            source=[];offset=[];total=[];centroids=[]
            for k in range(5):
                x=history['x_pixel'][k].ravel();y=history['y_pixel'][k].ravel()
                if k<4:
                    rec=recordmap[e,k];rec['response'].update(eta);q=targets[rec['name']];profile=np.exp(q[2])*rec['response'].values(x-q[0],y-q[1])
                    clone={**rec,'x':x,'y':y,'bounds':history['bounds_xyxy'][k].tolist()}
                    w=wcss[e]
                    def sky_at_detector(x,y,w=w):
                        ra,dec=w(x,y)[:2];return prior.tangent_plane_arcsec(np.stack([ra,dec],axis=-1),origin)
                    js.build_scene_projection(clone,sky_at_detector,anchors[k]);profile+=clone['scene_matrix']@sky[9*k:9*k+9]
                    center=np.mean([r['guess'] for r in records if r['exposure']==e],axis=0)
                    plane=np.c_[np.ones(len(x)),(x-center[0])/100,(y-center[1])/100]@backgrounds[e]
                    centroids.append(q[:2])
                else:
                    rec=fieldmap[e];rec['response'].update(eta);i=exposures.index(e);q=jointp[fieldstart+6*i:fieldstart+6*i+6]
                    profile=np.exp(q[2])*rec['response'].values(x-q[0],y-q[1]);o=rec['guess'];plane=np.c_[np.ones(len(x)),x-o[0],y-o[1]]@q[3:]
                    centroids.append(q[:2])
                source.append(profile.reshape(65,65));offset.append(plane.reshape(65,65));total.append((profile+plane).reshape(65,65))
            source=np.asarray(source);offset=np.asarray(offset);total=np.asarray(total);factor=history['flat']/float(history['photmjsr']);gain=history['gain_e_per_DN']
            path=OUT/f'{mode}_{e.replace("_cal.fits", "_pilot.npz")}'
            np.savez_compressed(path,labels=history['labels'],bounds_xyxy=history['bounds_xyxy'],x_pixel=history['x_pixel'],y_pixel=history['y_pixel'],centroids_native_xy=np.asarray(centroids),eta=eta,
                 positive_point_plus_coherent_sky_SCI=source,detector_plane_offset_SCI=offset,total_conditional_mean_SCI=total,
                 point_plus_sky_DN_s=source*factor,point_plus_sky_e_s=source*factor*gain,
                 total_mean_DN_s=total*factor,total_mean_e_s=total*factor*gain)
            mode_reports.append({'file':str(path.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'history_sha256':hashlib.sha256(historypath.read_bytes()).hexdigest(),'total_mean_min_SCI':float(total.min()),'offset_min_max_SCI':[float(offset.min()),float(offset.max())]})
        reports[mode]=mode_reports
    prior.dump(OUT/'receipt.json',{'scope':'Saved converged conditional first2GROUP opticalmean rendered on exactmatching65x65target+field history domains. No refit or unusedobservation. Trainingplusheld usesE1–3targettraining/all4field andindependentE4pointmeasurement.',
      'state_separation':'Positivepoint+coherentsky anddetectoroffset plane areseparate. A positive totalmean isnotproofthat thedetectoroffset isphotonillumination; consumers muststate orderive thatpartition beforecomputing λ-dependent noise orcollection. No negativeoffset isclipped into a photonrate.',
      'domain':'Finite25nativepixel positivePSF zerooutside; fourfinite.55arcsecsky patches/declaredGaussianbasis; largerhost/sourceandbackgroundphotonpartitionunresolved. ExactarchivedGWCS andsinglepixelintegratedPSFconvolution used.',
      'outputs':reports,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
    print('saved',sum(len(v) for v in reports.values()),'matching65pixelpilotmeans',flush=True)


if __name__=='__main__':main()

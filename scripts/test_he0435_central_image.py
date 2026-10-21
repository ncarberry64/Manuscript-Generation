"""Execute HST pixel tests, with no access to light-curve validation targets."""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"code"))
from strong_lensing_hst_central import (empirical_psf, fit_point, point_template,
    registered_light_center, robust_sigma, positivity_upper, nonnegative_profile,
    smooth_light_profile, assess_central_image_system, sky_offset_to_pixel)


def render_panels(panel,out):
    """Plot the saved numerical application, without changing any fits."""
    def draw(ax,a,title,residual,vmax,fontsize):
        if residual:
            im=ax.imshow(a,origin="lower",cmap="RdBu_r",vmin=-10,vmax=10)
        else:
            im=ax.imshow(a,origin="lower",norm=matplotlib.colors.AsinhNorm(.08,vmin=0,vmax=vmax))
        ax.set_title(title,fontsize=fontsize)
        ax.set_xticks([]);ax.set_yticks([])
        cb=plt.colorbar(im,ax=ax,fraction=.046,pad=.025)
        if residual:cb.set_ticks([-10,-5,0,5,10]);cb.set_label("residual / sigma",fontsize=fontsize-1)
        else:
            cb.set_ticks([v for v in [.01,.1,1,10,100] if v<=vmax])
            cb.set_label("electrons s$^{-1}$",fontsize=fontsize-1)
        cb.ax.tick_params(labelsize=fontsize-1)
    titles=["Observed central pixels","Positive galaxy + nucleus\n+ frozen central source",
            "Residual / marginal sigma","Predicted central light alone"]
    fig,axs=plt.subplots(2,4,figsize=(12,6),constrained_layout=True)
    for row,(band,d,m,r,c)in enumerate(panel):
        vmax=max(float(d.max()),float(c.max()))
        for col,a in enumerate([d,m,r,c]):draw(axs[row,col],a,f"{band}: {titles[col]}",col==2,vmax,8)
    fig.savefig(out/"hst_central_pixel_test.png",dpi=160)
    fig.savefig(out/"hst_central_pixel_test.pdf")
    plt.close(fig)
    band,d,m,r,c=next(p for p in panel if p[0]=="F814W")
    with plt.rc_context({"font.size":11}):
        fig,axs=plt.subplots(2,2,figsize=(7.2,6.2),constrained_layout=True)
        vmax=max(float(d.max()),float(c.max()))
        for k,a in enumerate([d,m,r,c]):draw(axs.flat[k],a,titles[k],k==2,vmax,11)
        fig.savefig(ROOT/"manuscript/figures/he0435_hst_central_pixel_test.pdf")
        fig.savefig(out/"hst_central_F814W_publication.png",dpi=180)
        plt.close(fig)


def render_saved_receipt():
    """Reapply saved fitted amplitudes for display; no optimization or refit."""
    out=ROOT/"artifacts/strong_lensing_revision_20261020"
    path=out/"hst_central_image.json";receipt=json.loads(path.read_text());panel=[]
    for band,old in receipt["bands"].items():
        scale=old["pixel_scale_arcsec"]
        with fits.open(ROOT/f"data/strong_lensing_revision_20261020/hst/{band}_lens_cutout.fits")as h:
            data=np.asarray(h["SCI"].data,float);w=np.asarray(h["WHT"].data,float)
        with fits.open(ROOT/f"data/strong_lensing_revision_20261020/hst/{band}_star2_cutout.fits")as h:
            psf,pc,_=empirical_psf(h["SCI"].data)
        sky=old["background_electrons_per_second"]
        err=np.sqrt(old["blank_sky_sigma_electrons_per_second"]**2*old["median_blank_WHT_exposure_like"]/np.maximum(w,1)+np.maximum(data-sky,0)/np.maximum(w,1))
        case=next(c for c in old["psf_cases"]if c["psf"]=="field_star_2")
        center=np.asarray(case["registered_light_center_pixel"]);cx,cy=np.rint(center).astype(int);half=int(np.ceil(.9/scale))
        sl=(slice(cy-half,cy+half+1),slice(cx-half,cx+half+1));patch=data[sl]-sky;sig=err[sl];local=center-[cx-half,cy-half]
        fit=next(f for f in old["smooth_and_nuclear_light_profile_fits"]if f["model"]=="flexible_cusp_plus_nucleus" and f["fixed_central_A_ratio"]>0)
        q=np.asarray(fit["parameters"]);yy,xx=np.indices(patch.shape);gx,gy=local+q[:2]/scale;co,si=np.cos(q[3]),np.sin(q[3])
        x=((xx-gx)*co+(yy-gy)*si)*scale;y=(-(xx-gx)*si+(yy-gy)*co)*scale
        kernel=point_template((25,25),[12,12],psf,pc);kernel/=kernel.sum()
        widths=np.unique(np.r_[.02,.08,.15,.3,.6,1.2])
        components=[fftconvolve(np.exp(-.5*((x/a)**2+(y/(a*q[2]))**2)),kernel,mode="same")for a in widths]
        components+=[point_template(patch.shape,[gx,gy],psf,pc),np.ones(patch.shape),(xx-local[0])*scale,(yy-local[1])*scale]
        pos=local+sky_offset_to_pixel(receipt["prediction"]["position_light_centered_arcsec"],scale,case["registration_rotation_rad"])
        central=point_template(patch.shape,pos,psf,pc)*case["fitted_A_flux_electrons_per_second"]*fit["fixed_central_A_ratio"]
        model=sum(a*c for a,c in zip(fit["amplitudes"],components))+central
        panel.append((band,patch,model,(patch-model)/sig,central))
    render_panels(panel,out)
    receipt["code_sha256"]["scripts/test_he0435_central_image.py"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    receipt["display_regeneration"]="Saved fitted parameters and amplitudes reapplied exactly; no science refit. Positive colour-bar zero tick removed; readable F814W publication panel added."
    path.write_text(json.dumps(receipt,indent=2)+"\n")


def assess_revised_predictions(input_path, output_path):
    """Pixel evaluation for a development candidate's full additional-image set."""
    incoming=json.loads(Path(input_path).read_text())
    original_path=ROOT/"artifacts/strong_lensing_revision_20261020/hst_central_image.json"
    original=json.loads(original_path.read_text())
    data_dir=ROOT/"data/strong_lensing_revision_20261020/hst"
    result={"status":"DEVELOPMENT_HST_CONSTRAINT_NOT_UNUSED_VALIDATION",
            "input_file":str(input_path),"input_sha256":hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
            "original_HST_receipt_sha256":hashlib.sha256(original_path.read_bytes()).hexdigest(),
            "photometric_condition":original["photometric_condition"],"models":{}}
    for model,spec in incoming["models"].items():
        images=spec["images"]
        positions=[i["position_light_centered_arcsec"]for i in images]
        ratios=[i["absolute_macro_flux_relative_A"]for i in images]
        by_band={}
        for band,old in original["bands"].items():
            scale=old["pixel_scale_arcsec"]
            with fits.open(data_dir/f"{band}_lens_cutout.fits")as h:
                data=np.asarray(h["SCI"].data,float);w=np.asarray(h["WHT"].data,float)
            sky=old["background_electrons_per_second"]
            variance=old["blank_sky_sigma_electrons_per_second"]**2*old["median_blank_WHT_exposure_like"]/np.maximum(w,1)+np.maximum(data-sky,0)/np.maximum(w,1)
            sigma=np.sqrt(variance)
            cases=[]
            for case in old["psf_cases"]:
                if case["psf"].startswith("field_star"):
                    n=case["psf"].split("_")[-1]
                    with fits.open(data_dir/f"{band}_star{n}_cutout.fits")as h:
                        psf,pc,_=empirical_psf(h["SCI"].data)
                else:
                    ax,ay=np.rint(case["quasar_fits"][0]["xy"]).astype(int)
                    half=int(np.ceil(.45/scale))
                    psf,pc,_=empirical_psf(data[ay-half:ay+half+1,ax-half:ax+half+1])
                def application(z,noise=1.,normal=.95):
                    assessment,_=assess_central_image_system(data-sky,noise*sigma,data.shape,
                        case["registered_light_center_pixel"],scale,positions,ratios,
                        case["fitted_A_flux_electrons_per_second"],psf,pc,z,normal,
                        rotation_rad=case["registration_rotation_rad"])
                    return assessment
                fixed=application([0,0]);best=fixed["upper_scale_95"]
                if best is None:
                    cases.append({"psf":case["psf"],"fixed_position":fixed,"free_centroid_upper_scale_95":None});continue
                zbest=[0,0]
                for dx in np.linspace(-.05,.05,11):
                    for dy in np.linspace(-.05,.05,11):
                        value=application([dx,dy])["upper_scale_95"]
                        if value>best:best=value;zbest=[dx,dy]
                opt=minimize(lambda z:-application(z)["upper_scale_95"],zbest,method="Powell",
                    bounds=[(-.05,.05),(-.05,.05)],options={"xtol":1e-7,"ftol":1e-8})
                if -opt.fun>best:best=float(-opt.fun);zbest=opt.x.tolist()
                stress=minimize(lambda z:-application(z,5.,.8)["upper_scale_95"],zbest,method="Powell",
                    bounds=[(-.05,.05),(-.05,.05)],options={"xtol":1e-7,"ftol":1e-8})
                cases.append({"psf":case["psf"],"fixed_position":fixed,
                    "free_centroid_upper_scale_95":best,
                    "fivefold_sigma_twenty_percent_lower_A_upper_scale_95":float(-stress.fun),
                    "centroid_offset_pixel_directions_arcsec":list(zbest)})
            valid=[c for c in cases if c["free_centroid_upper_scale_95"]is not None]
            by_band[band]={"cases":cases,"worst_psf_centroid_upper_scale_95":max((c["free_centroid_upper_scale_95"]for c in valid),default=None),
                "worst_adversarial_upper_scale_95":max((c["fivefold_sigma_twenty_percent_lower_A_upper_scale_95"]for c in valid),default=None)}
        result["models"][model]={"images":images,"number_additional_images":len(images),"bands":by_band}
    target=Path(output_path);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({m:{b:v["worst_psf_centroid_upper_scale_95"]for b,v in r["bands"].items()}for m,r in result["models"].items()},indent=2))


def main():
    data_dir = ROOT/"data/strong_lensing_revision_20261020/hst"
    out = ROOT/"artifacts/strong_lensing_revision_20261020"
    out.mkdir(parents=True, exist_ok=True)
    ast = json.loads((ROOT/"data/strong_lensing_observed_20261020/candidates/HE0435-1223/astrometry_redshifts.json").read_text())
    catalog = ast["astrometry"]["data"]
    g = np.array([catalog[4]["delta_alpha"], catalog[4]["delta_delta"]])
    sky = np.array([[r["delta_alpha"], r["delta_delta"]] for r in catalog[:4]])-g
    audit = json.loads((ROOT/"artifacts/strong_lensing_observed_20261020/macro_forward_audit.json").read_text())
    central = audit["conventional_census"]["images"][0]
    offset = np.array(central["position_light_centered_arcsec"])
    ratio = central["absolute_macro_flux_relative_A"]
    receipt = {"status":"EXECUTED_ORIGINAL_FROZEN_MACRO_CENTRAL_IMAGE_PIXEL_TEST",
               "prediction":central,"source_data_provenance":"data/strong_lensing_revision_20261020/hst/provenance.json",
               "conditional_R1_prediction":audit["conditional_R1_census"]["images"][0],
               "no_lightcurve_targets_consumed":True,
               "photometric_condition":"Observed central/A ratio equals frozen achromatic macro ratio; no differential extinction/microlensing/intrinsic-delay variability correction is derived here.",
               "bands":{},"code_sha256":{}}
    panel = []
    for band,scale,rough in [("F814W",0.05, np.array([80.,80.])),
                             ("F160W",0.12825000286,np.array([40.,40.]))]:
        with fits.open(data_dir/f"{band}_lens_cutout.fits") as h:
            data=np.asarray(h["SCI"].data,float); weight=np.asarray(h["WHT"].data,float)
            exp=float(h[0].header["EXPTIME"])
        yy,xx=np.indices(data.shape)
        rr=np.hypot(xx-rough[0],yy-rough[1])*scale
        blank=(rr>2.5)&(rr<3.7)&np.isfinite(data)&(weight>0.7*np.nanmedian(weight))
        sky_bg=float(np.median(data[blank])); sky_sigma=robust_sigma(data[blank])
        medweight=float(np.median(weight[blank]))
        variance=sky_sigma**2*medweight/np.maximum(weight,1)+np.maximum(data-sky_bg,0)/np.maximum(weight,1)
        sigma=np.sqrt(variance)
        cases=[]
        guesses=rough+np.c_[-sky[:,0],sky[:,1]]/scale
        # A coarse light-center guess differs slightly from the NIC2 reference.
        for i in range(4):
            gx,gy=np.rint(guesses[i]).astype(int)
            local=data[gy-3:gy+4,gx-3:gx+4]
            jy,jx=np.unravel_index(np.argmax(local),local.shape)
            guesses[i]=[gx-3+jx,gy-3+jy]
        templates=[]
        for star in [1,2,3]:
            with fits.open(data_dir/f"{band}_star{star}_cutout.fits") as h:
                psf,pc,bg=empirical_psf(h["SCI"].data)
            templates.append((f"field_star_{star}",psf,pc))
        # An A-derived template checks the different SED and detector sampling
        # of the very object whose flux normalizes the central prediction.
        ax,ay=np.rint(guesses[0]).astype(int)
        ahalf=int(np.ceil(.45/scale))
        apsf,apc,_=empirical_psf(data[ay-ahalf:ay+ahalf+1,ax-ahalf:ax+ahalf+1])
        templates.append(("quasar_A_empirical",apsf,apc))
        reference_images=None
        for label,psf,pc in templates:
            fits_q=[fit_point(data,sigma,q,psf,pc,0.30/scale) for q in guesses]
            center,rot,ast_resid=registered_light_center([q["xy"]for q in fits_q],sky,scale)
            a_flux=fits_q[0]["flux"]
            original_xy=center+sky_offset_to_pixel(offset,scale,rot)
            p=point_template(data.shape,original_xy,psf,pc)
            limit,meta=positivity_upper(data-sky_bg,sigma,p)
            # A conservative 5% aperture/PSF normalization nuisance acts against
            # exclusion. Its effect is reported separately from statistical noise.
            limit_ratio=limit/(0.95*a_flux)
            profile,nn,rrn=nonnegative_profile(data-sky_bg,sigma,p,ratio*a_flux)
            maxlim=0.;best_offset=None
            # Two-dimensional freedom of the lens-light centroid, wider than
            # the declared 3mas astrometric error, with no galaxy-centre prior.
            for dx in np.linspace(-.05,.05,11):
                for dy in np.linspace(-.05,.05,11):
                    pp=point_template(data.shape,original_xy+np.array([dx,dy])/scale,psf,pc)
                    lim,_=positivity_upper(data-sky_bg,sigma,pp)
                    if lim>maxlim:maxlim=lim;best_offset=[float(dx),float(dy)]
            def negative_limit(z, noise_factor=1):
                pp=point_template(data.shape,original_xy+z/scale,psf,pc)
                lim,_=positivity_upper(data-sky_bg,noise_factor*sigma,pp)
                return -lim
            continuous=minimize(negative_limit,best_offset,method="Powell",
                    bounds=[(-.05,.05),(-.05,.05)],options={"xtol":1e-7,"ftol":1e-8})
            if -continuous.fun>maxlim:
                maxlim=float(-continuous.fun);best_offset=continuous.x.tolist()
            stress=minimize(lambda z:negative_limit(z,5),best_offset,method="Powell",
                    bounds=[(-.05,.05),(-.05,.05)],options={"xtol":1e-7,"ftol":1e-8})
            case={"psf":label,"fitted_A_flux_electrons_per_second":a_flux,
                  "quasar_fits":[{k:(v.tolist()if isinstance(v,np.ndarray)else v)for k,v in q.items()}for q in fits_q],
                  "registered_light_center_pixel":center.tolist(),"registration_rotation_rad":rot,
                  "max_quasar_registration_residual_arcsec":ast_resid,
                  "central_pixel":original_xy.tolist(),"fixed_position_positive_light_95_upper_ratio":limit_ratio,
                  "free_centroid_50mas_square_95_upper_ratio":maxlim/(.95*a_flux),
                  "least_constraining_centroid_offset_pixel_directions_arcsec":best_offset,
                  "continuous_centroid_search":True,
                  "fivefold_sigma_twenty_percent_lower_A_normalization_upper_ratio":float(-stress.fun/(.8*a_flux)),
                  "upper_bound_details":meta,
                  "nonnegative_pixel_profile_chi2_at_frozen_ratio":profile,
                  "profile_chi2_is_not_a_detection_significance":True,
                  "minimum_required_extra_dimming_at_free_centroid_mag":max(0,float(-2.5*np.log10(maxlim/(.95*a_flux)/ratio)))}
            r1=audit["conditional_R1_census"]["images"][0]
            r1xy=center+sky_offset_to_pixel(r1["position_light_centered_arcsec"],scale,rot)
            r1p=point_template(data.shape,r1xy,psf,pc)
            r1lim,_=positivity_upper(data-sky_bg,sigma,r1p)
            case["conditional_R1_fixed_position_positive_light_95_upper_ratio"]=r1lim/(.95*a_flux)
            cases.append(case)
            if label=="field_star_2":reference_images=(data,sigma,psf,pc,center,a_flux,p,rot)
        data0,sigma0,psf0,pc0,center0,a_flux0,p0,rot0=reference_images
        # Central smooth-light nuisance fits are genuine applications. The
        # positivity test above remains valid for arbitrary stellar profiles,
        # host arcs, and any unresolved nonnegative nuclear component.
        half=int(np.ceil(.9/scale)); cx,cy=np.rint(center0).astype(int)
        sl=(slice(cy-half,cy+half+1),slice(cx-half,cx+half+1))
        patch=data0[sl]-sky_bg;err=sigma0[sl]; c_local=center0-[cx-half,cy-half]
        smooth=[]
        for min_width,nucleus,name in [(.08,False,"smooth_positive_MGE"),(.02,True,"flexible_cusp_plus_nucleus")]:
            previous=None
            for f in [0.,ratio]:
                fit,model=smooth_light_profile(patch,err,scale,psf0,pc0,c_local,offset,
                    f*a_flux0,minimum_width=min_width,nuclear_point=nucleus,start=previous,rotation_rad=rot0)
                previous=fit["parameters"]
                fit.update({"model":name,"fixed_central_A_ratio":f})
                smooth.append(fit)
                if name=="flexible_cusp_plus_nucleus" and f==ratio:
                    panel.append((band,patch,model,(patch-model)/err,p0[sl]*ratio*a_flux0))
        # Exact injection/recovery in the conservative nuisance class: its
        # profile upper limit shifts by the injected amplitude at fixed position.
        limit0,_=positivity_upper(data0-sky_bg,sigma0,p0)
        inject=[]
        rng=np.random.default_rng(4351223)
        core=p0>.03*p0.max()
        pcore=p0[core]
        # Deliberately worst-case boundary B=0 is a valid nuisance state and
        # tests marginal coverage of the conservative Bonferroni construction.
        sims=3000; covered=0
        for _ in range(sims):
            # Perfectly correlated Gaussian noise is a nontrivial check that
            # the bound does not require independent drizzled pixels.
            z=rng.normal(); d=ratio*a_flux0*pcore+z*sigma0[core]
            lim,_=positivity_upper(d,sigma0[core],pcore)
            covered+=lim>=ratio*a_flux0
        for f in [.01,ratio,.1]:
            injected=data0-sky_bg+f*a_flux0*p0
            lim,_=positivity_upper(injected,sigma0,p0)
            inject.append({"injected_ratio":f,"upper_bound_increment_ratio":(lim-limit0)/a_flux0,
                           "recovery_error_ratio":(lim-limit0)/a_flux0-f})
        receipt["bands"][band]={"pixel_scale_arcsec":scale,"exposure_seconds":exp,
             "background_electrons_per_second":sky_bg,"blank_sky_sigma_electrons_per_second":sky_sigma,
             "median_blank_WHT_exposure_like":medweight,
             "weight_formula":"variance = blank_sigma^2 median(WHT_blank)/WHT + max(SCI-sky,0)/WHT; WHT is exposure-like, never inverse variance.",
             "normalization_nuisance_fraction":.05,"centroid_search_square_half_width_arcsec":.05,
             "psf_cases":cases,"smooth_and_nuclear_light_profile_fits":smooth,
             "worst_psf_and_centroid_95_upper_ratio":max(c["free_centroid_50mas_square_95_upper_ratio"]for c in cases),
             "worst_psf_centroid_and_adversarial_noise_normalization_upper_ratio":max(c["fivefold_sigma_twenty_percent_lower_A_normalization_upper_ratio"]for c in cases),
             "injection_recovery":inject,
             "correlated_noise_boundary_simulation":{"seed":4351223,"simulations":sims,
                 "coverage":covered/sims,"noise":"perfectly correlated Gaussian pixel noise; B=0 boundary","nominal_coverage":.95}}
    render_panels(panel,out)
    for rel in ["code/strong_lensing_hst_central.py","scripts/test_he0435_central_image.py"]:
        receipt["code_sha256"][rel]=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
    (out/"hst_central_image.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({b:{"upper_ratio":v["worst_psf_and_centroid_95_upper_ratio"],
        "smooth_profiles":v["smooth_and_nuclear_light_profile_fits"]}for b,v in receipt["bands"].items()},indent=2))


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--prediction-json")
    parser.add_argument("--plot-only",action="store_true")
    parser.add_argument("--output",default=str(ROOT/"artifacts/strong_lensing_revision_20261020/hst_revised_central_image.json"))
    args=parser.parse_args()
    if args.plot_only:render_saved_receipt()
    elif args.prediction_json:assess_revised_predictions(args.prediction_json,args.output)
    else:main()

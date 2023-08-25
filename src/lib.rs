// Copyright 2023 Ahmad W. Huran
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use std::ffi::CStr;
use std::os::raw::{c_int, c_double};
use std::mem::MaybeUninit;

use anyhow::{Result, anyhow};

include!(concat!(env!("OUT_DIR"), "/libxc-bindings.rs"));



pub enum LibXCEvalType {
    EXC,
    VXC,
    FXC,
}

pub enum LibXCReturnType{
    EXC(Vec<f64>),

    LDAVXC(Vec<f64>),
    LDAFXC(Vec<f64>),

    GGAVXC(Vec<f64>, Vec<f64>),
    GGAFXC(Vec<f64>, Vec<f64>, Vec<f64>),

}

impl LibXCReturnType {
    pub fn get_exc(self) -> Option<Vec<f64>> {
        match self {
            Self::EXC(exc) => Some(exc),
            _ => None
        }
    }

    pub fn get_lda_vxc(self) -> Option<Vec<f64>> {
        match self {
            Self::LDAVXC(vxc) => Some(vxc),
            _ => None
        }
    }

    pub fn get_lda_fxc(self) -> Option<Vec<f64>> {
        match self {
            Self::LDAFXC(fxc) => Some(fxc),
            _ => None
        }
    }

    pub fn get_gga_vxc(self) -> Option<(Vec<f64>, Vec<f64>)> {
        match self {
            Self::GGAVXC(vrho, vsigma) => Some((vrho, vsigma)),
            _ => None
        }
    }

    pub fn get_gga_fxc(self) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        match self {
            Self::GGAFXC(v2rho2, v2rhosigma, v2sigma2) => Some((v2rho2, v2rhosigma, v2sigma2)),
            _ => None
        }
    }
}


#[derive(Debug)]
pub enum LibXCSpin{
    Unpolarized,
    Polarized,
}

#[derive(Debug)]
pub enum LibXCFamily{
    LDA,
    GGA,
    Other,
}


pub struct LibXCFunctional{
    id: u32,
    func: xc_func_type,
    spin: LibXCSpin,
    family: LibXCFamily,
    dim: Vec<usize>,

}

impl LibXCFunctional{

    pub fn new(id:u32, spin: LibXCSpin) -> Result<Self>{
        let func_id: c_int = id.try_into()?;
        let mut func: MaybeUninit<xc_func_type> = MaybeUninit::uninit();
        let _spin: c_int = match  spin {
            LibXCSpin::Unpolarized => XC_UNPOLARIZED.try_into()?,
            LibXCSpin::Polarized => XC_POLARIZED.try_into()?,

        };

        let err: c_int = unsafe {
            xc_func_init(func.as_mut_ptr(), func_id, _spin)
        };
        if err != 0 {
            return Err(
                anyhow!("Failed to initialize XC functiona id: {} and spin: {:?}",
                        id, spin
                        )
                )
        }
        let func = unsafe { func.assume_init()};

        let info: xc_func_info_type = unsafe{*func.info};
        let family_id = info.family as u32;
        let family = match family_id {
            XC_FAMILY_LDA => LibXCFamily::LDA,
            XC_FAMILY_GGA => LibXCFamily::GGA,
            _ =>  LibXCFamily::Other,

        };

        let dim: Vec<usize> = match (&spin, &family) {
            (LibXCSpin::Unpolarized, LibXCFamily::LDA) => vec![1, 1, 1],
            (LibXCSpin::Unpolarized, LibXCFamily::GGA) => vec![1, 1, 1],
            (LibXCSpin::Polarized, LibXCFamily::LDA) => vec![2, 2, 3],
            (LibXCSpin::Polarized, LibXCFamily::GGA) => vec![2, 3, 6],
            (_, LibXCFamily::Other) => vec![0, 0],
        };
        Ok(
            Self {
                id,
                func,
                spin,
                family,
                dim,
            }
        )
    }

    pub fn evaluate<S>(
        &self,
        rho: S,
        sigma: Option<S>,
        typ: LibXCEvalType,
        ) -> Result<LibXCReturnType>
    where S: AsRef<[f64]>
    {

        let rho_slice: &[f64] = rho.as_ref();
        let mut n_points = rho_slice.len();

        if n_points % self.dim[0] != 0 {
            return Err(anyhow!("Length of spin polarized rho must be  even!"))
        }

        n_points = n_points / self.dim[0];

        let func_ptr = &self.func as *const xc_func_type;
        let rho_ptr = rho_slice.as_ptr() as *const c_double;





        match (&self.family, &sigma) {
            (LibXCFamily::LDA, Some(_)) => return Err(
                anyhow!(
                    "Found 'sigma != None' while evaluating a \
                    'LibXCFamily::LDA' functional"
                    )
                ),

            (LibXCFamily::GGA, None) => return Err(
                anyhow!(
                    "Found 'sigma = None' while evaluating a \
                    'LibXCFamily::GGA' functional"
                    )
                ),


            (LibXCFamily::Other, _) => return Err(
                anyhow!(
                    "Only LDA and GGA are supported!"
                    )
                ),

            (_, _) => ()
        };

        match self.family {

            LibXCFamily::LDA => {
                match typ {
                    LibXCEvalType::EXC => {
                        let mut exc: Vec<f64> = vec![0.0; n_points];
                        let exc_ptr = exc.as_mut_ptr() as *mut c_double;
                        unsafe{
                            xc_lda_exc(
                                func_ptr,
                                n_points,
                                rho_ptr,
                                exc_ptr
                                )
                        };
                        return Ok(LibXCReturnType::EXC(exc))
                    },
                    LibXCEvalType::VXC => {
                        let mut vrho: Vec<f64> = vec![0.0; self.dim[1]*n_points];
                        let vrho_ptr = vrho.as_mut_ptr() as *mut c_double;
                        unsafe{
                            xc_lda_vxc(
                                func_ptr,
                                n_points,
                                rho_ptr,
                                vrho_ptr
                                )
                        };
                        return Ok(LibXCReturnType::LDAVXC(vrho))
                    },
                    LibXCEvalType::FXC => {
                        let mut v2rho2: Vec<f64> = vec![0.0; self.dim[2]*n_points];
                        let v2rho2_ptr = v2rho2.as_mut_ptr() as *mut c_double;
                        unsafe{
                            xc_lda_fxc(
                                func_ptr,
                                n_points,
                                rho_ptr,
                                v2rho2_ptr
                                )
                        };
                        return Ok(LibXCReturnType::LDAFXC(v2rho2))
                    },
                }
            },

            LibXCFamily::GGA => {
                let sig = sigma.unwrap();
                let sig_slice: &[f64] = sig.as_ref();

                if sig_slice.len() % self.dim[1] != 0 {
                    return Err(anyhow!("Length of spin polarized sigma must be 3 times the number of points!"))
                }
                let sig_ptr = sig_slice.as_ptr() as *const c_double;

                match typ {
                    LibXCEvalType::EXC => {
                        let mut exc: Vec<f64> = vec![0.0; n_points];
                        let exc_ptr = exc.as_mut_ptr() as *mut c_double;
                        unsafe{
                            xc_gga_exc(
                                func_ptr,
                                n_points,
                                rho_ptr,
                                sig_ptr,
                                exc_ptr
                                )
                        };
                        return Ok(LibXCReturnType::EXC(exc))
                    },
                    LibXCEvalType::VXC => {
                        let mut vrho: Vec<f64> = vec![0.0; self.dim[0]*n_points];
                        let vrho_ptr = vrho.as_mut_ptr() as *mut c_double;

                        let mut vsigma: Vec<f64> = vec![0.0; self.dim[1]*n_points];
                        let vsigma_ptr = vsigma.as_mut_ptr() as *mut c_double;

                        unsafe{
                            xc_gga_vxc(
                                func_ptr,
                                n_points,
                                rho_ptr,
                                sig_ptr,
                                vrho_ptr,
                                vsigma_ptr
                                )
                        };
                        return Ok(LibXCReturnType::GGAVXC(vrho, vsigma))
                    },
                    LibXCEvalType::FXC => {
                        let mut v2rho2: Vec<f64> = vec![0.0; self.dim[1]*n_points];
                        let v2rho2_ptr = v2rho2.as_mut_ptr() as *mut c_double;

                        let mut v2rhosigma: Vec<f64> = vec![0.0; self.dim[2]*n_points];
                        let v2rhosigma_ptr = v2rhosigma.as_mut_ptr() as *mut c_double;

                        let mut v2sigma2: Vec<f64> = vec![0.0; self.dim[2]*n_points];
                        let v2sigma2_ptr = v2sigma2.as_mut_ptr() as *mut c_double;

                        unsafe{xc_gga_fxc(
                                func_ptr,
                                n_points,
                                rho_ptr,
                                sig_ptr,
                                v2rho2_ptr,
                                v2rhosigma_ptr,
                                v2sigma2_ptr,
                                )
                        };
                        return Ok(
                            LibXCReturnType::GGAFXC(
                                v2rho2,
                                v2rhosigma,
                                v2sigma2)
                                  )
                    },
                };

            },

            _ => return Err(anyhow!("Only LDA and GGA are supported!"))
        }

    }

    pub fn get_family(&self) -> &LibXCFamily{
        return &self.family
    }

    pub fn get_id(&self) -> u32{
        return self.id
    }

    pub fn get_spin(&self) -> usize{
        match &self.spin {
            LibXCSpin::Unpolarized => 0,
            LibXCSpin::Polarized => 1,
        }
    }
}

impl Drop for LibXCFunctional{
    fn drop(&mut self) {
        let ptr = &mut self.func as *mut xc_func_type;
        unsafe{xc_func_end(ptr)}
    }
}

impl std::fmt::Display for LibXCFunctional {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let info: xc_func_info_type = unsafe{*self.func.info};
        let name = unsafe{CStr::from_ptr(info.name).to_str().unwrap()};
        let family_id = info.family as u32;
        let kind_id = info.kind as u32;

        let kind = match kind_id {
            XC_EXCHANGE => "an exchange functional",
            XC_CORRELATION => "a correlation functional",
            XC_EXCHANGE_CORRELATION => "an exchange-correlation functional",
            XC_KINETIC => "a kinetic energy functional",
            _ => "of unknown kind",
        };

        let family = match family_id {
            XC_FAMILY_LDA => "LDA",
            XC_FAMILY_GGA => "GGA",
            XC_FAMILY_HYB_GGA => "Hybrid GGA",
            XC_FAMILY_MGGA => "MGGA",
            XC_FAMILY_HYB_MGGA => "Hybrid MGGA",
            _ => "unknown"
        };

        let refs = info.refs;

        let mut ref_string = String::new();
        for (ii, ref_) in refs.clone().into_iter().enumerate(){
            let ix = ref_ as u32;
            if ix == 0{
                break;
            }
            let ref_ = unsafe{ *ref_};
            let ref_ = unsafe{ CStr::from_ptr(ref_.ref_).to_str().unwrap()};
            ref_string.push_str(format!("[{}] {}\n", ii+1, ref_).as_str());

        }




        write!(f, "The functional '{}' is {}, it belongs to the '{}' family and is defined in the reference(s):\n{}", name, kind, family, ref_string)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_lda_eval_exc_unpolarized() {
        let xx: Vec<f64> = vec![-1.0, 0.0, 1.0];
        let rho: Vec<f64> = xx.iter()
            .map(|xi| (-xi*xi).exp())
            .collect();
        let lda_c_pz = LibXCFunctional::new(9, LibXCSpin::Unpolarized).unwrap();
        let lda_x = LibXCFunctional::new(1, LibXCSpin::Unpolarized).unwrap();

        let e_c = lda_c_pz.evaluate(&rho, None, LibXCEvalType::EXC)
            .unwrap()
            .get_exc()
            .unwrap() ;
        let e_x = lda_x.evaluate(&rho, None, LibXCEvalType::EXC)
            .unwrap()
            .get_exc()
            .unwrap();

        let e_xc: Vec<f64> = e_c
            .iter()
            .zip(e_x)
            .map(|(c, x)| c+x)
            .collect();
        let ground_truth: Vec<f64> = vec![
            -0.5919756493441902,
            -0.8091965676851791,
            -0.5919756493441902,
        ];

        assert_eq!(e_xc, ground_truth)
    }

    #[test]
    fn test_lda_eval_exc_polarized() {
        let xx: Vec<f64> = vec![-1.0, 0.0, 1.0];
        let mut rho: Vec<f64> = Vec::with_capacity(6);
        for xi in xx {
            rho.push(0.5 * (-xi*xi).exp());
            rho.push(0.5 * (-xi*xi).exp());
        }


        let lda_c_pz = LibXCFunctional::new(9, LibXCSpin::Polarized).unwrap();
        let lda_x = LibXCFunctional::new(1, LibXCSpin::Polarized).unwrap();

        let e_c = lda_c_pz.evaluate(&rho, None, LibXCEvalType::EXC)
            .unwrap()
            .get_exc()
            .unwrap();

        let e_x = lda_x.evaluate(&rho, None, LibXCEvalType::EXC)
            .unwrap()
            .get_exc()
            .unwrap();

        let e_xc: Vec<f64> = e_c
            .iter()
            .zip(e_x)
            .map(|(c, x)| c+x)
            .collect();
        let ground_truth: Vec<f64> = vec![
            -0.5919756493441902,
            -0.8091965676851791,
            -0.5919756493441902,
        ];

        assert_eq!(e_xc, ground_truth)
    }

}

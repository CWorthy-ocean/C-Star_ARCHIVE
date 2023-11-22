%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINE_TRACERS: Return a list of tracers info to fill netcdf files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
    tracers= {...
	{'PO4'        ,'Phosphate'                                        ,'mMol PO4 m-3'   }
        {'NO3'        ,'Nitrate'                                          ,'mMol No3 m-3'   }
        {'SiO3'       ,'Silicate'                                         ,'mMol SiO3 m-3'  }
        {'NH4'        ,'Ammonia'                                          ,'mMol PO4 m-3'   }
        {'Fe'         ,'Iron'                                             ,'mMol Fe m-3'    }
        {'O2'         ,'Oxygen'                                           ,'mMol O2 m-3'    }
        {'DIC'        ,'Dissolved inorganic carbon'                       ,'mMol C m-3'     }
        {'Alk'        ,'Alkalinity'                                       ,'mMol m-3'       }
        {'DOC'        ,'Dissolved organic carbon'                         ,'mMol C m-3'     }
        {'DON'        ,'Dissolved organic nitrogen'                       ,'mMol N m-3'     }
        {'DOFE'       ,'Dissolved organic iron'                           ,'mMol Fe m-3'    }
        {'DOP'        ,'Dissolved organic phosphorus'                     ,'mMol P m-3'     }
        {'DOPR'       ,'refractory Dissolved organic phosphorus'          ,'mMol P m-3'     }
        {'DONR'       ,'refractory Dissolved organic nitrogen'            ,'mMol N m-3'     }
        {'ZOOC'       ,'Zooplankton carbon'                               ,'mg C m-3'       } % i.e. ZOOC=0.4*SPC=0.4*5*SPCHL
        {'SPC'        ,'Small phytoplankton carbon'                       ,'mMol C m-3'     } % i.e. 0.75*5*SPCHL (MF: initiate with lower Carbon!)
        {'SPCHL'      ,'Small phytoplankton chlorophyll'                  ,'mg Chla m-3'    } % for diaz_to_coocos: 0.45
        {'SPFE'       ,'Small phytoplankton iron'                         ,'mmol Fe m-3'    }
        {'SPCACO3'    ,'Small phytoplankton carbonate'                    ,'mmol CaCO3 m-3' } % i.e. 0.1*SPCHL
        {'DIATC'      ,'Diatom carbon'                                    ,'mMol C m-3'     } % i.e. DIATC=0.75*3*DIATCHL=0.75*3*0.1*SPCHL (MF: initiate with lower Carbon!)
        {'DIATCHL'    ,'Diatom chlorophyll'                               ,'mg Chla m-3'    } % i.e. 0.1*SPCHL=0.1*0.9*CHL
        {'DIATFE'     ,'Diatom iron'                                      ,'mmol Fe m-3'    } % i.e. 2e-5*DIATCHL
        {'DIATSI'     ,'Diatom silicate'                                  ,'mmol Si m-3 '   } % i.e. 0.1*DIATCHL
        {'DIAZC'      ,'Diazotroph carbon'                                ,'mmol C m-3'     } % for diaz_to_coocos: 5  (MF: initiate with lower Carbon!)
        {'DIAZCHL'    ,'Diazotroph chlorophyll'                           ,'mg Chla m-3'    } % i.e. DIAZCHL=0.01/0.9*SPCHL=0.01*CHL; for diaz_to_coocos: 1 (1*SPCHL)
        {'DIAZFE'     ,'Diazotroph iron'                                  ,'mmol Fe m-3'    } % for diaz_to_coocos: 1e-4*0.01/0.9
        {'NO2'        ,'Nitrite'                                          ,'mMol NO2 m-3'   }
        {'N2'         ,'Dinitrogen'                                       ,'mMol N2  m-3'   }
        {'N2O'        ,'Nitrous oxide'                                    ,'mMol N2O m-3'   }
        {'pCO2'       ,'Surface water pCO2'                               ,'uatm'           }
        {'CHLA'       ,'' ,'' }
        {'DIC_glodap' ,'' ,'' }
        {'Alk_glodap' ,'' ,'' }
        {'N2O_SIDEN'  ,'' ,'' }
        {'N2O_ATM'    ,'' ,'' }
        {'N2O_NEV'    ,'' ,'' }
        {'basindx'    ,'' ,'' }
        };

for trc=1:length(tracers)
    test = tracers{trc} ;
    bgctracers_list.name{trc}     = test{1} ;
    bgctracers_list.longname{trc} = test{2} ;
    bgctracers_list.units{trc}    = test{3} ;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINE_TRACERS: Return a list of tracers info to fill netcdf files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
    tracers= {...
	{'PO4'        ,'Dissolved Inorganic Phosphate'                                        ,'mMol PO4 m-3'   }
        {'NO3'        ,'Dissolved Inorganic Nitrate'                                          ,'mMol No3 m-3'   }
        {'SiO3'       ,'Dissolved Inorganic Silicate'                                         ,'mMol SiO3 m-3'  }
        {'NH4'        ,'Dissolved Ammonia'                                                    ,'mMol NH4 m-3'   }
        {'Fe'         ,'Dissolved Inorganic Iron'                                             ,'mMol Fe m-3'    }
        {'Lig'        ,'Iron Binding Ligand'                                                  ,'mMol Lig m-3'   }
        {'O2'         ,'Dissolved Oxygen'                                                     ,'mMol O2 m-3'    }
        {'DIC'        ,'Dissolved Inorganic Carbon'                                           ,'mMol C m-3'     }
        {'DIC_ALT_CO2','Dissolved Inorganic Carbon, Alternative CO2'                          ,'mMol C m-3'     }
        {'ALK'        ,'Alkalinity'                                                           ,'meq m-3'        }
        {'ALK_ALT_CO2','Alkalinity, Alternative CO2'                                          ,'meq m-3'        }
        {'DOC'        ,'Dissolved Organic Carbon'                                             ,'mMol C m-3'     }
        {'DON'        ,'Dissolved Organic Nitrogen'                                           ,'mMol N m-3'     }
        {'DOP'        ,'Dissolved Organic Phosphorus'                                         ,'mMol P m-3'     }
        {'DOPr'       ,'Refractory Dissolved Organic Phosphorus'                              ,'mMol P m-3'     }
        {'DONr'       ,'Refractory Dissolved Organic Nitrogen'                                ,'mMol N m-3'     }
        {'DOCr'       ,'Refractory Dissolved Organic Carbon'                                  ,'mMol C m-3'     }
        {'zooC'       ,'Zooplankton Carbon'                                                   ,'mg C m-3'       }
        {'spChl'      ,'Small Phytoplankton Chlorophyll'                                      ,'mg Chla m-3'    }
        {'spC'        ,'Small Phytoplankton Carbon'                                           ,'mMol C m-3'     }
        {'spP'        ,'Small Phytoplankton Phosphorous'                                      ,'mMol P m-3'     }
        {'spFe'       ,'Small Phytoplankton Iron'                                             ,'mmol Fe m-3'    }
        {'spCaCO3'    ,'Small Phytoplankton CaCO3'                                            ,'mmol CaCO3 m-3' }
        {'diatChl'    ,'Diatom Chlorophyll'                                                   ,'mg Chla m-3'    }
        {'diatC'      ,'Diatom Carbon'                                                        ,'mMol C m-3'     }
        {'diatP'      ,'Diatom Phosphorous'                                                   ,'mMol P m-3'     }
        {'diatFe'     ,'Diatom Iron'                                                          ,'mmol Fe m-3'    }
        {'diatSi'     ,'Diatom Silicate'                                                      ,'mmol Si m-3 '   }
        {'diazChl'    ,'Diazotroph Chlorophyll'                                               ,'mg Chla m-3'    }
        {'diazC'      ,'Diazotroph Carbon'                                                    ,'mmol C m-3'     }
        {'diazP'      ,'Diazotroph Phosphorous'                                               ,'mmol P m-3'     }
        {'diazFe'     ,'Diazotroph Iron'                                                      ,'mmol Fe m-3'    }
        };

for trc=1:length(tracers)
    test = tracers{trc} ;
    bgctracers_list.name{trc}     = test{1} ;
    bgctracers_list.longname{trc} = test{2} ;
    bgctracers_list.units{trc}    = test{3} ;
end





import cobra


def simulate_LB_media(model: cobra.Model) -> None:
    #lb medium ref: https://github.com/cdanielmachado/carveme/blob/master/carveme/data/input/media_db.tsv
    #TODO: what should the flux of all the metabolites here be set to?
    #TODO: simulated LB media here, but paper's media is '2 x YT'

    LB_MEDIA_COMP = ["EX_adn_e", "EX_ala__L_e", "EX_amp_e", "EX_arg__L_e",
                     "EX_aso3_e", "EX_asp__L_e", "EX_ca2_e", "EX_cbl1_e",
                     "EX_cd2_e", "EX_cl_e", "EX_cmp_e", "EX_cobalt2_e",
                     "EX_cro4_e", "EX_cu2_e", "EX_cys__L_e", "EX_dad_2_e",
                     "EX_dcyt_e", "EX_fe2_e", "EX_fe3_e", "EX_fol_e", "EX_glc__D_e",
                     "EX_glu__L_e", "EX_gly_e", "EX_gmp_e", "EX_gsn_e", "EX_h2o_e",
                     "EX_h2s_e", "EX_h_e", "EX_hg2_e", "EX_his__L_e", "EX_hxan_e",
                     "EX_ile__L_e", "EX_ins_e", "EX_k_e", "EX_leu__L_e", "EX_lipoate_e",
                     "EX_lys__L_e", "EX_met__L_e", "EX_mg2_e", "EX_mn2_e", "EX_mobd_e",
                     "EX_na1_e", "EX_nac_e", "EX_nh4_e", "EX_ni2_e", "EX_o2_e",
                     "EX_phe__L_e", "EX_pheme_e", "EX_pi_e", "EX_pnto__R_e",
                     "EX_pro__L_e", "EX_pydx_e", "EX_ribflv_e", "EX_ser__L_e",
                     "EX_so4_e", "EX_thm_e", "EX_thr__L_e", "EX_thymd_e", "EX_trp__L_e",
                     "EX_tyr__L_e", "EX_ump_e", "EX_ura_e", "EX_uri_e", "EX_val__L_e",
                     "EX_zn2_e"
                      ]

    for metabolite in LB_MEDIA_COMP:
        model.medium[metabolite] = 10



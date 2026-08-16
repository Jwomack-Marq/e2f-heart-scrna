function state(; kwargs...)
	defaults = (
		CCND = 0, # Cyclin D: Cyclin D involved in G1 to S transition.
		CDK4 = 1, # Cyclin-Dependent Kinase 4: Cyclin-dependent kinase associated with CCND.
		CCND_CDK4 = 0, # Active CCND-CDK4 Complex.
		CCND_CDK4p = 0, # Inactive CCND-CDK4 Complex.
		E2F = 0, # Transcription factor activating genes for cell cycle progression.
		RB1 = 1, # Retinoblastoma Protein: Tumor suppressor protein regulating cell cycle entry.
		RB_E2F = 0, # Retinoblastoma Protein-E2F Complex: Complex inhibiting E2F activity.
		pRB_E2F = 0, # Phosphorylated Retinoblastoma Protein-E2F Complex: Partially activated E2F.
		ppRB = 0, # Hyper-phosphorylated Retinoblastoma Protein: Phosphorylated form released from E2F.
		CCNE = 0, # Cyclin E: Cyclin E involved in G1 to S transition.
		CDK2 = 2, # Cyclin-Dependent Kinase 2: Cyclin-dependent kinase associated with CCNE and CCNA.
		CCNE_CDK2 = 0, # Active CCNE-CDK2 Complex.
		CCNE_CDK2p = 0, # Inactive CCNE-CDK2 Complex.
		CDC25A = 0, # Inactive, dephosphorylated form of phosphatase Cdc25A.
		CDC25Ap = 0, # Active, phosphorylated form of phosphatase Cdc25A.
		CCNA = 0, # Cyclin A: Cyclin A involved in S phase and G2 progression.
		CCNA_CDK2 = 0, # Active CCNA-CDK2 Complex.
		CCNA_CDK2p = 0, # Inactive CCNA-CDK2 Complex.
		CCNB = 0, # Cyclin B1: Regulatory protein involved in mitosis, forming a complex with CDK1 to become active and initiate mitosis.
		CDK1 = 1, # Cyclin-dependent kinase 1 regulates the cell cycle by phosphorylating target proteins.
		CCNB_CDK1 = 0, # Maturation Promoting Factor, also known as MPF, is a complex of cyclin B1 and CDK1 that drives cells into mitosis.
		CCNB_CDK1p = 0, # Pre-Maturation Promoting Factor, the inactive form of MPF.
		WEE1 = 0.5, # Kinase WEE1 phosphorylates CDK1, inhibiting its activity and delaying entry into mitosis.
		WEE1p = 0.5, # Phosphorylated form of WEE1 (inactive).
		CDC25C = 1, # Inactive CDC25C phosphatase.
		CDC25Cp = 0, # Phosphorylated form of CDC25C, active in promoting entry into mitosis.
		p21 = 0, # p21, also known as cyclin-dependent kinase inhibitor 1, regulates cell cycle progression by inhibiting CDK activity.
		p21_CCNB_CDK1 = 0, # Complex of p21 with MPF (Cyclin B1/CDK1), regulating entry into mitosis.
		PLK1 = 0, # Polo-like kinase 1, involved in various stages of mitosis, including centrosome maturation and spindle assembly (inactive).
		PLK1p = 0, # Phosphorylated form of PLK1, regulating its activity (active).
		APCC = 0.5, # Anaphase-Promoting Complex/Cyclosome without specific activator.
		APCCP = 0.5, # Anaphase-Promoting Complex/Cyclosome activated by CDC20.
		PPase = 1, # Protein phosphatase involved in cell cycle regulation.
		PPaseP = 0, # Phosphorylated form of PPase, inactive.
		APCC_CDH1 = 0, # Anaphase-Promoting Complex/Cyclosome activated by CDH1.
		APCCP_CDH1 = 0, # inactive version of APCC_CDH1 that has been phosphorylated.
		CDH1 = 1, # Subunit of APC/C, activating it during late mitosis and G1 phase.
		CDH1p = 0, # Phosphorylated form of CDH1, inactivating it.
		CDC20 = 0, # Activator of APC/C, initiating anaphase by targeting proteins for degradation.
		CDC20p = 0, # Phosphorylated form of CDC20, possibly regulating its function.
		APCCP_CDC20 = 0, # Anaphase-Promoting Complex/Cyclosome activated by CDC20, targeting proteins for degradation in mitosis.
		PTTG1 = 1, # Pituitary Tumor Transforming Gene 1, involved in cell cycle regulation and tumorigenesis.
		PTTG1p = 0, # Phosphorylated form of PTTG1, possibly modulating its activity.
		LMNA = 1, # Lamin A/C, a structural protein of the nuclear envelope involved in maintaining nuclear shape and rigidity.
		LMNAp = 0, # Phosphorylated form of LMNA, possibly affecting its role in nuclear structure.
		EMI1 = 0, # EMI1, a negative regulator of the APC/C, preventing its activation.
		EMI1p = 0, # Phosphorylated form of EMI1 (inactive).
		APCCP_CDC20_EMI1 = 0, # APC/C activated by CDC20 and inhibited by EMI1, preventing degradation of specific substrates.
		APCC_CDH1_EMI1 = 0, # APC/C activated by CDH1 and inhibited by EMI1, preventing degradation of specific substrates.
		SCF = 0, # Skp, Cullin, F-box containing Complex: Regulates cell cycle proteins through ubiquitin.
		p21_CCNE_CDK2 = 0, # Inhibited CCNE_CDK2:CDK4 complex bonded to p21.
		p21_CCNA_CDK2 = 0, # Inhibited CCNA:CDK2 complex bonded to p21.
		ATM = 1, # Ataxia Telangiectasia Mutated, a kinase involved in DNA damage response.
		ATMp = 0, # Phosphorylated form of ATM, active in DNA damage response.
		Chk2 = 1, # Checkpoint kinase 2, involved in DNA damage response.
		Chk2p = 0, # Phosphorylated form of Chk2, active in DNA damage response.
		p53 = 0, # Protein complex that activates DNA damage pathway (inactive).
		p53p = 0, # Protein complex that activates DNA damage pathway (active).
		MDM2 = 0, # E3 ubiquitin-protein ligase that regulates p53 activity.
		p53_MDM2 = 0, # Complex of p53 and MDM2, regulating p53 activity.
		CDT1 = 0, # NA
		Geminin = 0, # NA
		Geminin_CDT1 = 0 # NA
		)	
		user_state = merge(defaults, kwargs)
        return ComponentVector{Float64}(user_state)
end

function state_names()
	[
		"CCND",
		"CDK4",
		"CCND_CDK4",
		"CCND_CDK4p",
		"E2F",
		"RB1",
		"RB_E2F",
		"pRB_E2F",
		"ppRB",
		"CCNE",
		"CDK2",
		"CCNE_CDK2",
		"CCNE_CDK2p",
		"CDC25A",
		"CDC25Ap",
		"CCNA",
		"CCNA_CDK2",
		"CCNA_CDK2p",
		"CCNB",
		"CDK1",
		"CCNB_CDK1",
		"CCNB_CDK1p",
		"WEE1",
		"WEE1p",
		"CDC25C",
		"CDC25Cp",
		"p21",
		"p21_CCNB_CDK1",
		"PLK1",
		"PLK1p",
		"APCC",
		"APCCP",
		"PPase",
		"PPaseP",
		"APCC_CDH1",
		"APCCP_CDH1",
		"CDH1",
		"CDH1p",
		"CDC20",
		"CDC20p",
		"APCCP_CDC20",
		"PTTG1",
		"PTTG1p",
		"LMNA",
		"LMNAp",
		"EMI1",
		"EMI1p",
		"APCCP_CDC20_EMI1",
		"APCC_CDH1_EMI1",
		"SCF", 
		"p21_CCNE_CDK2",
		"p21_CCNA_CDK2",
		"ATM",
		"ATMp",
		"Chk2", 
		"Chk2p" ,
		"p53", 
		"p53p",
		"MDM2",
		"p53_MDM2",
		"CDT1",
		"Geminin",
		"Geminin_CDT1"
	]
end
workflow PRSmix {

	Int disk = 100
	Int ncores = 2
	Int memory = 8

    call s1_harmonize_SNPeffects {
        input: 
            disk = disk,
            ncores = ncores,
            memory = memory
    }

    call s2_computePRS {
        input: 
            weight_file = s1_harmonize_SNPeffects.weight_out_file,

            disk = disk,
            ncores = ncores,
            memory = memory
    }


    call s3_combine_PRS {
        input: 
            original_beta_files_list = s1_harmonize_SNPeffects.weight_out_file,
            score_files_list = s2_computePRS.score_out,

            disk = disk,
            ncores = ncores,
            memory = memory
    }


    meta {
        author : "Buu Truong"
        email : "btruong@broadinstitute.org"
        description : "This workflow runs PRSmix - see the README on the github for more information - https://github.com/buutrg/PRSmix"
    }

	output {
        File weight_out_file = s1_harmonize_SNPeffects.weight_out_file
        File score_file = s2_computePRS.score_out
        File prsmix_output = s3_combine_PRS.prsmix_output
	}

}

task s1_harmonize_SNPeffects {

    Boolean a_runStep1 = false
    File? s1_weight_file

	File? ref_file
	File? pgs_folder_tar 
	File? pgs_list
	Int snp_col 
	Int a1_col 
	Int beta_col
	Boolean isheader = true
	Int chunk_size
	String out = "harmonized_snpeff.txt"

	Int ncores
	Int disk
	Int memory


	command <<<

        # echo "ref_file:"$ref_file
        # if [ "${ref_file}"=="" ]; then ref_file="-1"; fi
        # if [ "${pgs_folder_tar}"=="" ]; then pgs_folder_tar="-1"; fi
        # if [ "${pgs_list}"=="" ]; then pgs_list="-1"; fi

        # echo "ref_file:"$ref_file
        # echo "args: ${ref_file} ${pgs_folder_tar} ${pgs_list} ${ncores} ${snp_col} ${a1_col} ${beta_col} ${isheader} ${chunk_size} ${out} ${a_runStep1} ${s1_weight_file}"

        # echo ${a_runStep1}

        R --no-save --args ${ref_file} ${pgs_folder_tar} ${pgs_list} ${ncores} ${snp_col} ${a1_col} ${beta_col} ${isheader} ${chunk_size} ${out} ${a_runStep1} ${s1_weight_file} << RSCRIPT

            library(PRSmix)
            library(data.table)
            options(datatable.fread.datatable=FALSE)

            print("s1_harmonize_SNPeffects")

            args = commandArgs(TRUE)
            ls()
            print(args)

            ref_file = args[1]
            pgs_folder_tar = args[2]
            pgs_list = args[3]
            ncores = as.numeric(args[4])
            snp_col = as.numeric(args[5])
            a1_col = as.numeric(args[6])
            beta_col = as.numeric(args[7])
            isheader = as.logical(args[8])
            chunk_size = as.numeric(args[9])
            out = args[10]
            a_runStep1 = as.logical(args[11])
            s1_weight_file = args[12]

            print(a_runStep1)
            if (!a_runStep1) {
                system(paste0("cp ", s1_weight_file, " ", out))
                q()
            }

            ls()
            print(ref_file)
            print(pgs_folder_tar)

            system(paste0("tar -xvf ", pgs_folder_tar))
            system("ls")
            system("ls score_folder")
            pgs_folder = "score_folder"

            system(paste0("ls ", pgs_folder))

            harmonize_snpeffect_toALT(
                ref_file = ref_file, 
                pgs_folder = pgs_folder,
                pgs_list = pgs_list,
                snp_col = snp_col,
                a1_col = a1_col,
                beta_col = beta_col,
                isheader = isheader,
                chunk_size = chunk_size,
                ncores = ncores,
                out = out 
            )

        RSCRIPT

    >>>

	output {
		File weight_out_file = out
	}


	runtime {
		docker: "buutrg/prsmix:1.0.0"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${ncores}"
		bootDiskSizeGb: 100
	}

}

task s2_computePRS {

    Boolean a_runStep2 = false
    File? score_inp

	File bed 
	File bim 
	File fam 
	File weight_file
	String out

	Int disk = 30
	Int ncores = 1
	Int memory = 30

	command <<<

		R --no-save --args ${bed} ${weight_file} ${out} ${a_runStep2} ${score_inp} << RSCRIPT


			library(PRSmix)
            library(data.table)
            options(datatable.fread.datatable=FALSE)

			print("s2_computePRS")
            args = commandArgs(TRUE)
            ls()
            print(args)

            geno = substring(args[1], 1, nchar(args[1]) - 4)
            weight_file = args[2]
            out = args[3]
            a_runStep2 = as.logical(args[4])
            score_inp = args[5]

            print(a_runStep2)
            if (!a_runStep2) {
                system(paste0("cp ", score_inp, " ", out, ".sscore"))
                q()
            }

			compute_PRS(geno = geno, weight_file = weight_file, out = out, plink2_path="/usr/bin/plink2")

		RSCRIPT

    >>>

	output {
		File score_out = "${out}.sscore"
	}

	runtime {
		docker: "buutrg/prsmix:1.0.0"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${ncores}"
		bootDiskSizeGb: 100
	}

}

task s3_combine_PRS {

	File pheno_file
    File covariate_file
	File score_files_list
	File trait_specific_score_file
	String pheno_name
	Boolean isbinary
	String out
	Boolean liabilityR2
	String IID_pheno
	String covar_list
	String cat_covar_list = "none"
	Boolean is_extract_adjSNPeff = false
	File original_beta_files_list
	String train_size_list
	String power_thres_list
	String pval_thres_list
	Boolean read_pred_training
	Boolean read_pred_testing

	Int ncores
    Int disk = 30
	Int memory = 30


	command <<<

        R --no-save --args ${pheno_file} ${covariate_file} ${score_files_list} ${trait_specific_score_file} ${pheno_name} ${out} ${isbinary} ${liabilityR2} ${IID_pheno} ${covar_list} ${cat_covar_list} ${ncores} ${is_extract_adjSNPeff} ${original_beta_files_list} ${train_size_list} ${power_thres_list} ${pval_thres_list} ${read_pred_training} ${read_pred_testing} << RSCRIPT

			library(PRSmix)
            library(data.table)
            options(datatable.fread.datatable=FALSE)

			print("combine PRS")
            args = commandArgs(TRUE)
            ls()
            print(args)

            allPGS_list = NULL
            metascore = NULL
            liabilityR2 = F
            IID_pheno = "IID"
            covar_list = c("age", "sex", paste0("PC", 1:10))
            cat_covar_list = NULL
            ncores = 1
            is_extract_adjSNPeff = F
            original_beta_files_list = NULL
            train_size_list = NULL
            training_result_file = NULL
            power_thres_list = c(0.95)
            pval_thres_list = c(0.05)
            nfold_cv = 3
            read_pred_training = FALSE
            read_pred_testing = FALSE
            debug = F

            # args = c("pheno.txt", "covar.txt", "allscores.sscore", "trait_specific_score_file.txt", "pheno", "cPRS", "false", "false", "IID", "covar1,covar2,covar3,covar4,covar5,covar6,covar7,covar8,covar9,covar10,covar11,covar12", "none", "2", "false", "map_dir.txt", "300", "0.95", "0.05", "false", "false")

            pheno_file = args[1]
            covariate_file = args[2]
            score_files_list = args[3]
            trait_specific_score_file = args[4]
            pheno_name = args[5]
            out = args[6]
            isbinary = as.logical(args[7])
            liabilityR2 = as.logical(args[8])
            IID_pheno = args[9]
            covar_list = args[10]
            cat_covar_list = args[11]
            ncores = as.numeric(args[12])
            is_extract_adjSNPeff = as.logical(args[13])
            original_beta_files_list = args[14]
            train_size_list = args[15]
            power_thres_list = args[16]
            pval_thres_list = args[17]
            read_pred_training = as.logical(args[18])
            read_pred_testing = as.logical(args[19])

            if (cat_covar_list == "none") { 
                cat_covar_list = NULL
            } else {
                cat_covar_list = unlist(strsplit(cat_covar_list, split=","))
            }

            covar_list = unlist(strsplit(covar_list, split=","))
            pval_thres_list = as.numeric(unlist(strsplit(pval_thres_list, split=",")))
            power_thres_list = as.numeric(unlist(strsplit(power_thres_list, split=",")))
            train_size_list = as.numeric(unlist(strsplit(train_size_list, split=",")))

            combine_PRS(
                pheno_file = pheno_file,
                covariate_file = covariate_file,
                score_files_list = score_files_list,
                trait_specific_score_file = trait_specific_score_file,
                pheno_name = pheno_name,
                out = out,
                isbinary = isbinary,
                liabilityR2 = liabilityR2,
                IID_pheno = IID_pheno,
                covar_list = covar_list,
                cat_covar_list = cat_covar_list,
                ncores = ncores,
                is_extract_adjSNPeff = is_extract_adjSNPeff,
                original_beta_files_list = original_beta_files_list,
                train_size_list = train_size_list,
                power_thres_list = power_thres_list,
                pval_thres_list = pval_thres_list,
                read_pred_training = read_pred_training, 
                read_pred_testing = read_pred_testing
                )

            system("ls")
            system("mkdir prsmix_output")

            system("mv *_train_allPRS.txt prsmix_output")
            system("mv *_test_allPRS.txt prsmix_output")

            system("mv *_train_df.txt prsmix_output")
            system("mv *_test_df.txt prsmix_output")

            system("mv *_test_summary_traitPRS_withPRSmix.txt prsmix_output")
            system("mv *_test_summary_traitPRS_withPRSmixPlus.txt prsmix_output")

            system("mv *_time_PRSmix.txt prsmix_output")
            system("mv *_time_PRSmixPlus.txt prsmix_output")

            system("tar -cvf prsmix_output.tar prsmix_output")

		RSCRIPT
	>>>

	output {
        File prsmix_output = "prsmix_output.tar"
	}

	runtime {
		docker: "buutrg/prsmix:1.0.0"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${ncores}"
		bootDiskSizeGb: 50
    }

}

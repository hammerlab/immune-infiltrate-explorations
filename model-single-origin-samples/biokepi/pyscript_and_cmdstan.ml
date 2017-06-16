#use "topfind";;
#thread

#use "/tmp/biokepi_machine.ml";;

open Biokepi
open KEDSL
open Cmdliner

(* simple aliases *)
let host = Ketrew.EDSL.Host.parse "/tmp/KT-coclomachine/"
let run_program = Machine.run_program biokepi_machine
let (//) = Filename.concat
let install_path = install_tools_path // "conda"

let ii_dir = "/nfs-pool/ii"
let ii_conda_requirements_path = ii_dir // "conda_requirements.txt"
let pip_requirements_path = ii_dir // "pip_requirements.txt"
let reference_build = "b38"

let ii_home_dir = "/modelcache/eliza-immune/immune-infiltrate-explorations/model-single-origin-samples/"
(* let stan_model_path = ii_home_dir // stan_model *)
let model_dir = ii_home_dir // "models"
let rdump_dir = ii_home_dir // "rdump-data"
let output_dir = ii_home_dir // "model_output"

let tmp_data_file_witness = rdump_dir // "model6.2.data.R"

(* Create a new Conda env *)
let conda_env =
  Setup.Conda.(
    setup_environment
      ~python_version:`Python3
      install_path
      "ii"
  )
let conda_ensure = Setup.Conda.(configured ~run_program ~host ~conda_env)
let conda_init = Setup.Conda.init_env ~conda_env ()
(* end of Conda setup *)

(*node that installs all immune infiltrate tools into ii environment, and
  creates empty file to check if installation exists*)
let ii_witness = install_path // ".fullii"
let ii_node =
  workflow_node (single_file ~host ii_witness)
    ~name:"Install immune infiltrate conda environment"
    ~edges: [
      depends_on conda_ensure;
    ]
    ~make:(run_program
             Program.(
               conda_init
               && sh "ls -al /nfs-pool"
               && shf "ls -al %s" ii_dir
               && shf "ls -al %s" ii_conda_requirements_path
               && shf "conda install --file=%s" ii_conda_requirements_path
               && shf "pip install -r %s" pip_requirements_path
               && shf "touch %s" ii_witness
             )
          )
(* end of ii install *)

(* Runs the python script that should create {modelname}.data.R as witness file *)
let python_rdata_node
    python_script
    =
  workflow_node (single_file ~host tmp_data_file_witness)
    ~name:("Execute python script, " ^ (Filename.basename python_script))
    ~edges:[
      depends_on (Biokepi.Tools.Pyensembl.cache_genome ~run_with:biokepi_machine ~reference_build);
      depends_on ii_node;
      depends_on conda_ensure;
    ]
    ~make:(run_program
             Program.(
               conda_init (*initialize conda environment*)
               && Biokepi.Tools.Pyensembl.(set_cache_dir_command ~run_with:biokepi_machine)
               && shf "python %s" python_script
             )
          )

(* Build and train the model *)

let submit_job
  python_script
  stan_model
  =
  let master_node =
    workflow_node without_product
    ~name:("Run rdump script " ^ python_script ^ ", train model " ^ stan_model)
    ~edges:[
       depends_on (python_rdata_node python_script);
        depends_on (
          Biokepi.Tools.Cmdstan.(fit_model
                                              ~stan_model:(model_dir // stan_model)
                                              ~fit_method:"variational"
					      ~data_file:tmp_data_file_witness
                                              (* ~data_file:(rdump_dir // (stan_model ^ ".data.R")) *)
                                              ~output_file:(output_dir // (stan_model ^ "_output.csv"))
                                              ~run_with:biokepi_machine));
      ]
  in
  Ketrew.Client.submit_workflow master_node

(*Command line options*)
let stan_model =
  let doc = "Name of stan model to build within immune-infiltrate model directory, without .stan suffix" in
  Arg.(required & pos 1 (some string) None & info [] ~docv:"STAN_MODEL" ~doc)

let python_script =
  let doc = "Path to .py you want to execute" in
  Arg.(required & pos 0 (some string) None & info [] ~docv:"PYTHON_FILE" ~doc)

let cmd =
  let doc = "Extract and save rdump; train stan model on it" in
  let version = "0.0.0" in
  let man = [
    `S "Description";
    `P "$(tname) gets rdata in python and trains stan model via cmdstan";
  ] in
  Term.(const
          submit_job
        $ python_script
        $ stan_model
       ),
  Term.(info "get rdata, train model" ~version ~doc ~man)

let () =
  match Cmdliner.Term.eval cmd with
  | `Error _ -> exit 1
  | _ -> exit 0

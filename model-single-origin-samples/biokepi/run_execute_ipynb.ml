#use "topfind";;
#thread

#use "biokepi_machine.ml";;

open Biokepi
open KEDSL
open Cmdliner

(* simple aliases *)
let host = Ketrew.EDSL.Host.parse "/tmp/KT-coclomachine/"
let run_program = Machine.run_program biokepi_machine
let (//) = Filename.concat
let install_path = install_tools_path // "conda"
let ii_conda_requirements_path = "/nfs-pool/ii/conda_requirements.txt"
let pip_requirements_path = "/nfs-pool/ii/pip_requirements.txt"
let reference_build = "b38"

(* Create a new Conda env *)
let conda_env =
  Setup.Conda.(
    setup_environment
      ~python_version:`Python3
      install_path
      "fullii"
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
               && shf "conda install --file=%s" ii_conda_requirements_path
               && shf "pip install -r %s" pip_requirements_path
               && shf "touch %s" ii_witness
             )
          )
(* end of ii install *)

(*Submit the task*)
let submit_job
    my_notebook
    timeout =
  let master_node =
    workflow_node without_product
      ~name:("Execute jupyter notebook, " ^ (Filename.basename my_notebook))
      ~edges:[
        depends_on (Biokepi.Tools.Pyensembl.cache_genome ~run_with:biokepi_machine ~reference_build);
        depends_on ii_node;
        depends_on conda_ensure;
      ]
      ~make:(run_program
               Program.(
                 conda_init (*initialize conda environment*)
                 && Biokepi.Tools.Pyensembl.(set_cache_dir_command ~run_with:biokepi_machine)
                 && shf "jupyter nbconvert --to notebook --execute %s --ExecutePreprocessor.timeout=%d" my_notebook timeout
               )
            )
  in
  Ketrew.Client.submit_workflow master_node

(*Command line options*)
let my_notebook =
  let doc = "Path to Jupyter notebook you want to execute" in
  Arg.(required & pos 0 (some string) None & info [] ~docv:"JUPYTER_IPYNB" ~doc)

let timeout =
  let doc="The maximum time (in seconds) each notebook cell is allowed to run. If the execution takes longer an exception will be raised. Can set to -1 for no maximum timeout." in
  Arg.(required & pos 1 (some int) None & info [] ~docv:"TIMEOUT_SECS" ~doc)

let cmd =
  let doc = "Execute jupyter notebook. Saves to {mynotebook}.nbconvert.ipynb" in
  let version = "0.0.0" in
  let man = [
    `S "Description";
    `P "$(tname) runs jupyter nbconvert.";
  ] in
  Term.(const
    submit_job
        $ my_notebook
        $ timeout
       ),
  Term.(info "run_nbconvert" ~version ~doc ~man)

let () =
  match Cmdliner.Term.eval cmd with
  | `Error _ -> exit 1
  | _ -> exit 0

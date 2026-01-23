

using JuMP
import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF


function save_model(model, instance_name)
    mof = FF.Model(format = FF.FORMAT_MOF)
    MOI.copy_to(mof, JuMP.backend(model)) 
    MOI.write_to_file(mof, "$instance_name.mof.json")
end


function load_model(instance_name)
    mof = FF.Model(format = FF.FORMAT_MOF)
    MOI.read_from_file(mof, "$instance_name.mof.json")
end

function build_qp_from_mof(model)
    
end


#
# speed comparison with JuMP AD
#

function report(m)

    eval_f = 0.0
    eval_g = 0.0
    eval_grad_f = 0.0
    eval_jac_g = 0.0
    eval_hesslag = 0.0

    jeval_f = 0.0
    jeval_g = 0.0
    jeval_grad_f = 0.0
    jeval_jac_g = 0.0
    jeval_hesslag = 0.0

    TapeInterface.reset_timer(m.internalModel.evaluator)
    for i=1:5
    	solve(m)
    	eval_f += m.internalModel.evaluator.eval_f_timer
        eval_g += m.internalModel.evaluator.eval_g_timer
        eval_grad_f += m.internalModel.evaluator.eval_grad_f_timer
        eval_jac_g += m.internalModel.evaluator.eval_jac_g_timer
        eval_hesslag += m.internalModel.evaluator.eval_hesslag_timer
       
       	jeval_f += m.internalModel.evaluator.jeval_f_timer
        jeval_g += m.internalModel.evaluator.jeval_g_timer
        jeval_grad_f += m.internalModel.evaluator.jeval_grad_f_timer
        jeval_jac_g += m.internalModel.evaluator.jeval_jac_g_timer
        jeval_hesslag += m.internalModel.evaluator.jeval_hesslag_timer
    end


    ratio_f = eval_f/jeval_f
    ratio_g = eval_g/jeval_g
    ratio_grad_f = eval_grad_f/jeval_grad_f
    ratio_jac_g = eval_jac_g/jeval_jac_g
    ratio_hesslag = eval_hesslag/jeval_hesslag


    # @show ratio_f2
    # @show ratio_g2
    # @show ratio_grad_f2
    # @show ratio_jac_g2
    # @show ratio_hesslag2

return  ((eval_f,  eval_g,  eval_grad_f,  eval_jac_g,  eval_hesslag),
        (jeval_f, jeval_g, jeval_grad_f, jeval_jac_g, jeval_hesslag),
        (ratio_f, ratio_g, ratio_grad_f, ratio_jac_g, ratio_hesslag))

# jeval_f = m.internalModel.evaluator.jd.eval_f_timer
# jeval_g = m.internalModel.evaluator.jd.eval_g_timer
# jeval_grad_f = m.internalModel.evaluator.jd.eval_grad_f_timer
# jeval_jac_g = m.internalModel.evaluator.jd.eval_jac_g_timer
# jeval_hesslag = m.internalModel.evaluator.jd.eval_hesslag_timer
# ratio_f = eval_f/jeval_f
# ratio_g = eval_g/jeval_g
# ratio_grad_f = eval_grad_f/jeval_grad_f
# ratio_jac_g = eval_jac_g/jeval_jac_g
# ratio_hesslag = eval_hesslag/jeval_hesslag
# @show ratio_f
# @show ratio_g
# @show ratio_grad_f
# @show ratio_jac_g
# @show ratio_hesslag

end

# TIMES = 2

matrix = Array{Any,2}(6,4)
matrix[1,1:4] =["         ","ReverseDiffTape", "JuMP internal", "ratio" ]
matrix[2,1] = "eval_f" 
matrix[3,1] = "eval_g" 
matrix[4,1] = "eval_grad_f" 
matrix[5,1] = "eval_jac_g"
matrix[6,1] = "eval_hesslag"  

# for i=2:1:TIMES
(t,jt,r) = report(m)
for i in 1:length(t)
    matrix[i+1,2] = t[i]
end
for i in 1:length(jt)
    matrix[i+1,3] = jt[i]
end
for i in 1:length(t)
    matrix[i+1,4] = r[i]
end

# end
@show matrix



struct AMGCLInfo
    iters::Cint
    residual::Cdouble
end

myjson(param)=JSON3.write(param)
myjson(::Nothing)=""
myjson(s::String)=s

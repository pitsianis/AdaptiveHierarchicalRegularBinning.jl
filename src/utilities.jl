function enable_multithreading()
  Core.eval(@__MODULE__, :(multithreading() = true))
end

function disable_multithreading()
  Core.eval(@__MODULE__, :(multithreading() = false))
end

debugmode() = false

function enable_debugging()
  Core.eval(@__MODULE__, :(debugmode() = true))
end

function disable_debugging()
  Core.eval(@__MODULE__, :(debugmode() = false))
end

macro threading(test)
    esc(:(if $(@__MODULE__).multithreading()
      Threads.@threads($test)
    else
      $test
    end))
end

macro debugging(test)
    esc(:(if $(@__MODULE__).debugmode()
      $test
    end))
end

macro autoinfiltrate(cond=true)
  esc(:(if $(@__MODULE__).debugmode()
    if isdefined(Main, :Infiltrator) && $cond
      Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
    end
  end))
end
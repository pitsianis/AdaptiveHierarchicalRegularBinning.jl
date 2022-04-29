begin
  multithreading() = true

  macro threading(test)
      esc(:(if $(@__MODULE__).multithreading()
        Threads.@threads($test)
      else
        $test
      end))
  end
end

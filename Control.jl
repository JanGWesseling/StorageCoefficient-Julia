module Control

  include("Types.jl")
  include("DataManager.jl")
  include("VanGenuchten.jl")

  using DataFrames
  using Printf
  using Random
  using Plots
  using StatsBase
  using Statistics
  using Dates

  staring = Array{Main.Control.Types.MVGParam}
  bofekProfile = Array{Main.Control.Types.BofekProfile}
  flux = Array{Main.Control.Types.Flux}
  sample = Array{Main.Control.Types.Sample}(undef,1)
  combination = Array{Int64}(undef, 1, 1)
  soilLayer = Array{Main.Control.Types.SoilLayer}(undef,1)
  storageVolume = Array{Float64}(undef,1,1)
  zForStaring = Array{Float64}(undef,111)
  vForStaring = Array{Float64}(undef,111)
  swapOutput = Array{Main.Control.Types.SwapOutput2}(undef,1)
  critical = Array{Main.Control.Types.Critical}(undef,1)
  newCritical = Array{Main.Control.Types.Critical}(undef,1)

  swapDir = "/home/wesseling/DataDisk/Wesseling/WesW/Projects/Achilles/Swap/"
  drz = 10.0
  test = false

  function createArrayWithParams(aLayers :: Array{Main.Control.Types.SoilLayer}, aPeriod :: String)
    try
      try
        global combination = Array{Int64}(undef,1,1)
        nLayers = size(aLayers,1)
        resize!(sample, nLayers)
        samples = Array{Int64}(undef, nLayers)
        for i in 1:nLayers
          myMVG = Main.Control.Types.MVGParam(-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
          myArray = Array{Main.Control.Types.MVGParam}(undef,1)
          myArray[1] = myMVG
          mySample = Main.Control.Types.Sample(myArray)
          mySample.param = DataManager.readSampleData(aLayers[i].staringId, aPeriod)
          sample[i] = mySample
          samples[i] = size(sample[i].param,1)
        end
        nLayers = size(samples,1)
        nSamples = 1
        for i in 1:nLayers
          nSamples *= samples[i]
       end

       if nSamples > 0
         global combination = Array{Int64}(undef, nSamples, nLayers)
         for i in 1:nLayers
           nSame = nSamples
           for l in 1:i
             nSame ÷= samples[l]
           end
           nTimes = 1
           if i > 1
             for l in 1:i-1
               nTimes *= samples[l]
             end
           end
           nRepeat = nSamples ÷ (nSame * nTimes)
           pos = 0
           for l in 1:nTimes
             for j in 1:nRepeat
               for k in 1:nSame
                 pos += 1
                 global combination[pos,i] = j
               end
             end
           end
         end
       end
      catch ex
        println("????ERROR in createArrayWithParams: ", ex)
      end
    finally
    end
  end

  function storageForInfiltration(aGwl :: Float64, aFlux :: Main.Control.Types.Flux)
    storage = 0.0
    try
      try
        layer = 1
        while aGwl < soilLayer[layer].bottom && layer < size(soilLayer,1)
          layer += 1
        end

        counter = 1
        maxHeight = abs(aGwl) - drz
        z = 1.0;
        dH = -0.1
        dHOld = 0.0
        h = -1.0
        dZ = 0.010
        k = VanGenuchten.conductivity(soilLayer[layer].param, -1.0)
        theta = VanGenuchten.moistureContent(soilLayer[layer].param, -0.5)
        storage = (1.0 + 0.5 * dZ) * (soilLayer[layer].param.thetaS - theta)
        while z < maxHeight && k > 1.0e-16 && layer >= 1  && h > -1.0e8
          dHOld = -1.0e6
          dH = 0.0
          j = 0
          while j < 1000 && abs(dH - dHOld) > 1.0e-6
            j += 1;
            dHOld = dH;
            k = VanGenuchten.conductivity(soilLayer[layer].param, h - 0.5 * (dH + dHOld))
            if k > 1.0e-12
              dH = -dZ * (1.0 + ( aFlux.value / k ))
            end
          end
          h += 0.5 * (dH + dHOld)
          if h > -1.0
            h = -1.0
          end
          z = z + dZ
          theta = VanGenuchten.moistureContent(soilLayer[layer].param, h)
          storage += dZ * (soilLayer[layer].param.thetaS - theta)

          counter += 1

          if aGwl + z > soilLayer[layer].top
            layer -= 1
          end
        end
      catch ex
        println("???ERROR in storageForInfiltration: ",ex)
      end
    finally
    end
    return storage
  end

  function storageForRise(aGwl :: Float64, aFlux :: Main.Control.Types.Flux)
    storage = 0.0
    try
      try
        layer = 1
        while aGwl < soilLayer[layer].bottom && layer < size(soilLayer,1)
          layer += 1
        end

        maxHeight = abs(aGwl) - drz
        h = -1.0
        z = 1.0
        dh = -0.1
        dz = 0.10
        k = 1.0
        hold = 0.0
        zold = 0.0
        thetaSat = 0.0
        topReached = false
        counter = 0
        while k > 1.0e-16 && dz > 1.0e-5 && layer >= 0 && h > -1.0e8 && (!topReached)
          hold = h
          zold = z
          if h < -1000.0
            dh = -1.0
          else
            if h < -10000.0
              dh = -10.0
            else
              if h < -100000.0
                dh = -100.0
              end
            end
          end
          h = h + dh
          k = VanGenuchten.conductivity(soilLayer[layer].param, h - 0.5 * dh)
          if k > 1.0e-12
            dz = -dh / (1.0 + (aFlux.value / k))
          end
          z += dz
          thetaSat = soilLayer[layer].param.thetaS

          nextLayer = false
          if aGwl + z >= soilLayer[layer].top
            slope = dz/dh
            h = h - dh + (soilLayer[layer].top - (aGwl + z + dz)) / slope
            z = soilLayer[layer].top - aGwl
            nextLayer = true
          end
          theta = VanGenuchten.moistureContent(soilLayer[layer].param, 0.5 * (h + hold))
          storage += (z - zold) * (thetaSat - theta)
          counter += 1

          if nextLayer
            layer -= 1
            nextLayer = false
          end
          if z >= maxHeight
            topReached = true;
          end
        end
        if !topReached
          storage += thetaSat * (soilLayer[layer].top - aGwl - z)
          if layer > 1
            for i in 1:layer-1
              storage += (soilLayer[i].param.thetaS - soilLayer[i].param.thetaR) * soilLayer[i].size
            end
          end
        end
      catch ex
        println("???ERROR in storageForRise: ",ex)
      end
    finally
    end
    return storage
  end

  function assignParameterValues(aPos :: Int64)
    try
      try
        for i in 1:size(soilLayer,1)
          global soilLayer[i].param = sample[i].param[combination[aPos,i]]
        end
      catch ex
        println("???ERROR in assignParamValues: ", ex)
      end
    finally
    end
  end


  function generateValues(aFlux :: Main.Control.Types.Flux, aNumber :: Int64, aPeriod :: String)
    try
      try
        nCombi = size(combination,1)
        if aPeriod == "new"
          global storageVolume = Array{Float64}(undef,12,nCombi)
          for i in 1:nCombi
            assignParameterValues(i)
            for j in 1:12
              gwl = -10.0 * j
              stor = 0.0
              if aFlux.value < 0.0
                stor = 10.0 * storageForInfiltration(gwl,aFlux)
              else
                stor = 10.0 * storageForRise(gwl,aFlux)
              end
              global storageVolume[j,i] = stor
            end
          end
        else
          global storageVolume = Array{Float64}(undef,12,aNumber)
          seed = round(Int64, 1000*time())
          rng = MersenneTwister(seed)
          for i in 1:aNumber
            if mod(i,1000) == 0
              println(i)
            end
            pos = rand(rng, 1:nCombi)
            assignParameterValues(pos)
            for j in 1:12
              gwl = -10.0 * j
              stor = 0.0
              if aFlux.value < 0.0
                stor = 10.0 * storageForInfiltration(gwl,aFlux)
              else
                stor = 10.0 * storageForRise(gwl,aFlux)
              end
              global storageVolume[j,i] = stor
            end
          end
        end
      catch ex
        println("???ERROR in generateValues: ", ex)
      end
    finally
    end
  end

  function computeWithStaring(aFlux :: Main.Control.Types.Flux)
    try
      try
        for i in 1:size(soilLayer,1)
          soilLayer[i].param = staring[soilLayer[i].staringId]
        end
        for i in 1:111
          gwl = -1.0 * (i + 9)
          stor = 0.0
          if aFlux.value < 0.0
            stor = 10.0 * storageForInfiltration(gwl,aFlux)
          else
            stor = 10.0 * storageForRise(gwl,aFlux)
          end
#          println(gwl,stor)
          global zForStaring[i] = -1.0 * abs(gwl)
          global vForStaring[i] = stor
        end
      catch ex
        println("???ERROR in computeWithStaring: ", ex)
      end
    finally
    end
  end

  function plotStorage(aBofekId :: Int64, aFlux :: Main.Control.Types.Flux, aNumber :: Int64, aPeriod :: String)
    try
      try
        nValues = size(storageVolume,2)
        myTitle = "Profiel=" * string(aBofekId) * ", q=" * aFlux.name * ", n=" * string(nValues)
        depth = 0
        nVals = size(storageVolume,1) * size(storageVolume,2)
        val = Array{Float64}(undef,nVals)
        x = Array{Int64}(undef,nVals)
        n = 0
        ave = Array{Float64}(undef,size(storageVolume,1))
        med = Array{Float64}(undef,size(storageVolume,1))
        q10 = Array{Float64}(undef,size(storageVolume,1))
        q90 = Array{Float64}(undef,size(storageVolume,1))
        dep = Array{Float64}(undef,size(storageVolume,1))
        v = Array{Float64}(undef,size(storageVolume,2))
        for i in 1:size(storageVolume,1)
          depth = -10 * i
          for j in 1:size(storageVolume,2)
            n += 1
            x[n] = depth
            val[n] = storageVolume[i,j]
            v[j] = storageVolume[i,j]
          end
          dep[i] = depth
          ave[i] = mean(v)
          med[i] = median(v)
          q = quantile(v,[0.1,0.9])
          q10[i] = q[1]
          q90[i] = q[2]
        end
        p = plot(xmirror=true, legend=:topright, ylab="Grondwaterdiepte (cm)", xlab="Berging (mm)", title=myTitle)
        p = plot!(p,val,x,seriestype=:scatter, markercolor=:blue, linewidth=0, markershape=:circle, markersize=5, label="Berekend")
        p = plot!(p,q90,dep,seriestype=:scatter, markercolor=:lime, linewidth=0, markershape=:square, markersize=6, label = "10%")
        p = plot!(p,q10,dep,seriestype=:scatter, markercolor=:aqua, linewidth=0, markershape=:square, markersize=6, label = "90%")
        p = plot!(p,ave,dep,seriestype=:scatter, markercolor=:orange, linewidth=0, markershape=:octagon, markersize=8, label = "Gemiddelde")
        p = plot!(p,med,dep,seriestype=:scatter, markercolor=:yellow, linewidth=0, markershape=:rtriangle, markersize=8, label = "Mediaan")

        p = plot!(p, vForStaring, zForStaring, seriestype=:line, line=:solid, color=:red, label="Staring")
        display(p)
        fileName = "/home/wesseling/DataDisk/Wesseling/WesW/Projects/Achilles/Plots/storage_" * string(aBofekId) * "_" * string(aFlux.id) * "_" * string(aNumber) * "_" * aPeriod * ".png"
        savefig(p, fileName)
      catch ex
        println("???ERROR in plotStorage", ex)
      end
    finally
    end
  end

  function createSoilLayers(aLocation)
    try
      try
        if aLocation == "nba"
          # create soil for Nieuw-Beertha
          resize!(soilLayer,5)
          myMVG = Main.Control.Types.MVGParam(-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
          soilLayer[1] = Main.Control.Types.SoilLayer(1,8,0,-8,11,myMVG)
          soilLayer[2] = Main.Control.Types.SoilLayer(2,12,-8,-20,11,myMVG)
          soilLayer[3] = Main.Control.Types.SoilLayer(3,30,-20,-50,30,myMVG)
          soilLayer[4] = Main.Control.Types.SoilLayer(4,20,-50,-70,35,myMVG)
          soilLayer[5] = Main.Control.Types.SoilLayer(5,430,-70,-500,35,myMVG)
        else
          if aLocation == "kat"
            resize!(soilLayer,5)
            myMVG = Main.Control.Types.MVGParam(-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
            soilLayer[1] = Main.Control.Types.SoilLayer(1,10,0,-10,9,myMVG)
            soilLayer[2] = Main.Control.Types.SoilLayer(2,15,-10,-25,9,myMVG)
            soilLayer[3] = Main.Control.Types.SoilLayer(3,25,-25,-50,28,myMVG)
            soilLayer[4] = Main.Control.Types.SoilLayer(4,20,-50,-70,28,myMVG)
            soilLayer[5] = Main.Control.Types.SoilLayer(5,230,-70,-300,28,myMVG)
          else
            # create B2/O2 soil
            resize!(soilLayer,5)
            myMVG = Main.Control.Types.MVGParam(-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
            soilLayer[1] = Main.Control.Types.SoilLayer(1,8,0,-30,2,myMVG)
            soilLayer[2] = Main.Control.Types.SoilLayer(2,12.0,-30,-400,20,myMVG)
          end
        end
      catch ex
        println("???ERROR in createSoilLayers: ",ex)
      end
    finally
    end
  end

  function createSoilPhysicsString(aRun :: Int64, aPos :: Int64)
    soils = ""
    mySample = [9,9,28,28,28]
    try
      try
        for i in 1:size(soilLayer,1)
          if aRun < 0
            if aRun > -10
              global soilLayer[i].param = staring[soilLayer[i].staringId]
            else
              sampleData = DataManager.readSampleData2(mySample[i])
              global soilLayer[i].param = sampleData[1]
            end
          else
            global soilLayer[i].param = sample[i].param[combination[aPos,i]]
          end
          l = soilLayer[i].param.l
          if l < -25.0
            l = -24.99
          end
          if l > 25.0
            l = 24.99
          end
          kSat = soilLayer[i].param.kSat
          if kSat < 1.0
            kSat += 1.0
          end
          soils *= string(i) * "  "
          soils *= @sprintf("%.3f", soilLayer[i].param.thetaR) * "  "
          soils *= @sprintf("%.3f", soilLayer[i].param.thetaS) * "  "
          soils *= @sprintf("%.6f", soilLayer[i].param.alpha) * "  "
          soils *= @sprintf("%.6f", soilLayer[i].param.n) * "  "
          soils *= @sprintf("%.3f", kSat) * "  "
          soils *= @sprintf("%.3f", l) * "  "
          soils *= @sprintf("%.6f", soilLayer[i].param.alpha) * "  0.0   "
          soils *= @sprintf("%.3f", kSat) * "  "
          soils *= "1315.0 \n"
        end
      catch ex
        println("???ERROR in createSoilPhysics: ",ex)
      end
    finally
    end
    return soils
  end

  function performSwapRun(aRun :: Int64, aSoil :: Int64, aLocation :: String)
    try
      try
        df = DateFormat("dd-u-yyyy HH:MM:SS")
       	startTime = "  started at " * Dates.format(Dates.now(), df)
        println(" run ", aRun, "  soil ", aSoil, startTime)
        soilString = createSoilPhysicsString(aRun, aSoil)
        Swap.createSwpFile(swapDir, soilString, aLocation)
        Swap.runSwap(swapDir, aLocation)
        if aLocation == "nba" || aLocation == "kat"
          out = Swap.getOutputData2(swapDir, aLocation)
          DataManager.storeSwapData(aRun, out, aLocation)
       else
         global swapOutput[aRun] = Swap.getOutputData(swapDir)
       end
      catch ex
        println("???ERROR in performSwapRun: ", ex)
      end
    finally
    end
  end

  function testWithSwap(aLocation)
    try
      try
        draw = false
        createSoilLayers(aLocation)
        createArrayWithParams(soilLayer, "new")
        println(size(combination,1), " possible combinations")
        if draw
#          DataManager.deleteFromNBA()
          n1 = -99
          n2 = -99
          nCombi = size(combination,1)
          seed = round(Int64, 1000*time())
          rng = MersenneTwister(seed)
          for i in n1:n2
            pos = rand(rng, 1:nCombi)
            performSwapRun(i, pos, aLocation)
          end
        else
          n = size(combination,1)
          n = 3
          resize!(swapOutput,n)
          for i in 1:n
            println(i)
            performSwapRun(i,i,aLocation)
          end
        end
        # run with Staring units
#        run = -1
#        performSwapRun(run,run,aLocation)
      catch ex
        println("???ERROR in testWithSwap: ",ex)
      end
    finally
    end
  end

  function exploreKsat()
    out = "id & min & max & ave & median & 5 & 95 & staring \n"
    try
      try
        for i in 1:size(staring,1)
          myMVG = Main.Control.Types.MVGParam(-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
          myArray = Array{Main.Control.Types.MVGParam}(undef,1)
          myArray[1] = myMVG
          mySample = Main.Control.Types.Sample(myArray)
          mySample.param = DataManager.readSampleData(staring[i].id,"all")
          sample[1] = mySample
          x = Array{Float64}(undef, size(sample[1].param,1))
          if size(sample[1].param,1) > 0
            for j in 1:size(sample[1].param,1)
              x[j] = sample[1].param[j].kSat
            end
            ave = mean(x)
            geo = geomean(x)
            med = median(x)
            q = quantile(x,[0.05,0.95])
            mi = minimum(x)
            ma = maximum(x)
            out *= string(staring[i].id) * " & "
            out *= @sprintf("%.2f", mi) * " & "
            out *= @sprintf("%.2f", ma) * " & "
            out *= @sprintf("%.2f", ave) * " & "
            out *= @sprintf("%.2f", geo) * " & "
            out *= @sprintf("%.2f", med) * " & "
            out *= @sprintf("%.2f", q[1]) * " & "
            out *= @sprintf("%.2f", q[2]) * " & "
            out *= @sprintf("%.2f", staring[i].kSat) * "\n"
          else
            out *= string(staring[i].id) * " & - & - & - & - & - & - & - & " * @sprintf("%.2f", staring[i].kSat) * "\n"
          end
        end
        outfil = open("/home/wesseling/DataDisk/Wesseling/WesW/Projects/Achilles/Doc/ksat.txt","w")
        write(outfil,out)
        close(outfil)
        println(out)

      catch ex
          println("???ERROR in exploreKsat: ", ex)
      end
    finally
    end
  end

  function generateZValues(aFlux :: Float64, aNumber :: Int64, aPeriod :: String)
    try
      try
        if aPeriod == "new"
          n = size(combination,1)
          if n > 0
            if test && n > aNumber
              n = 100
            end
            resize!(newCritical, n)
            for i in 1:n
              if mod(i,100) == 0
                println(i)
              end
              assignParameterValues(i)
              global newCritical[i] = computeCriticalValues(aFlux)
            end
          end
        else
          resize!(critical,aNumber)
          nCombi = size(combination,1)
          seed = round(Int64, 1000*time())
          rng = MersenneTwister(seed)
          for i in 1:aNumber
            if mod(i,1000) == 0
              println(i)
            end
            pos = rand(rng, 1:nCombi)
            assignParameterValues(pos)
            global critical[i] = computeCriticalValues(aFlux)
          end
        end
      catch ex
        println("???ERROR in generateZValues: ", ex)
      end
    finally
    end
  end

  function computeCriticalValues(aFlux :: Float64)
    critValue = Main.Control.Types.Critical(0.0,0.0)
    storage = 0.0
    try
      try
        drz = -10.0
        volumeRootZone = 0.0
        layer = 1
        while drz < soilLayer[layer].bottom
          theta = VanGenuchten.moistureContent(soilLayer[layer].param, -16000.0)
          thetaSat = soilLayer[layer].param.thetaS
          volumeRootZone += soilLayer[layer].size * (thetaSat - theta)
          layer += 1
        end
        theta = VanGenuchten.moistureContent(soilLayer[layer].param, -16000.0)
        thetaSat = soilLayer[layer].param.thetaS
        if layer > 1
          volumeRootZone += (soilLayer[layer-1].bottom - drz) * (thetaSat - theta)
        else
          volumeRootZone -= drz * (thetaSat - theta)
        end

        gwl = -20.0
        nextGwl = true
        gwl2 = -20.0
        gwl1 = soilLayer[size(soilLayer,1)].bottom
        while nextGwl
          gwl = round(0.5 * (gwl1 + gwl2))
          storage = 0.0
#          println("gwl1=",gwl1," gwl2=",gwl2," gwl=",gwl)
          layer = 1
          while gwl < soilLayer[layer].bottom && layer < size(soilLayer,1)
            layer += 1
          end

          maxHeight = drz - gwl
          h = -1.0
          z = 1.0
          dh = -0.1
          dz = 0.10
          k = 1.0
          hold = 0.0
          zold = 0.0
          thetaSat = 0.0
          topReached = false
          counter = 0
          while k > 1.0e-16 && dz > 1.0e-5 && layer >= 0 && h > -16000 && (!topReached)
            hold = h
            zold = z
            if h < -1000.0
              dh = -1.0
            else
              if h < -10000.0
                dh = -10.0
              else
                if h < -100000.0
                  dh = -100.0
                end
              end
            end
            h = h + dh
            k = VanGenuchten.conductivity(soilLayer[layer].param, h - 0.5 * dh)
            if k > 1.0e-12
              dz = -dh / (1.0 + (aFlux / k))
            end
            z += dz
            thetaSat = soilLayer[layer].param.thetaS

            nextLayer = false
            if gwl + z >= soilLayer[layer].top
              slope = dz/dh
              h = h - dh + (soilLayer[layer].top - (gwl + z + dz)) / slope
              z = soilLayer[layer].top - gwl
              nextLayer = true
            end
            theta = VanGenuchten.moistureContent(soilLayer[layer].param, 0.5 * (h + hold))
            storage += (z - zold) * (thetaSat - theta)
            counter += 1

            if nextLayer
              layer -= 1
              nextLayer = false
            end
            if z >= maxHeight
              topReached = true;
            end
          end
          if topReached
            gwl2 = gwl
          else
            gwl1 = gwl
          end
          if abs(gwl2 - gwl1) < 1.1
            nextGwl = false
            critValue.distance = gwl
            critValue.volume = 10.0 * (storage + volumeRootZone)
          end
        end
      catch ex
        println("???ERROR in computeCritValues: ",ex)
      end
    finally
    end
    return critValue
  end

  function plotCriticalValues(aCritical :: Main.Control.Types.Critical, aBofekId :: Int64, aFlux :: Float64, aNumber :: Int64)
    try
      try
        m = 0
        if size(newCritical,1) > 1
          m = size(newCritical,1)
        end
        myTitle = "Profiel=" * string(aBofekId) * ", q=" * @sprintf("%.1f",10*aFlux) * " mm/d, n=" * string(aNumber) * ", m=" * string(m)
        p = plot(legend=:topright, xlab="Kritieke grondwaterstand (cm)", ylab="Kritieke berging (mm)", title=myTitle,width=1200,height=800,xlim=(-200.0,0.0), ylim=(0.0,300.0))
        x = Array{Float64}(undef,1)
        y = Array{Float64}(undef,1)
        nVals = size(critical,1)
        if nVals > 1
          resize!(x,nVals)
          resize!(y,nVals)
          for i in 1:nVals
            x[i] = critical[i].distance
            y[i] = critical[i].volume
          end
          p = plot!(p,x,y,seriestype=:scatter, markercolor=:blue, linewidth=0, markershape=:circle, markersize=5, label="Oud")
        end

        if size(newCritical,1) > 1
          nVals = size(newCritical,1)
          resize!(x,nVals)
          resize!(y,nVals)
          for i in 1:nVals
            x[i] = newCritical[i].distance
            y[i] = newCritical[i].volume
          end
          p = plot!(p,x,y,seriestype=:scatter, markercolor=:lime, linewidth=0, markershape=:diamond, markersize=5, label="Nieuw")
        end

        resize!(x,1)
        resize!(y,1)
        x[1] = aCritical.distance
        y[1] = aCritical.volume
        p = plot!(p,x,y,seriestype=:scatter, markercolor=:red, linewidth=0, markershape=:square, markersize=6, label = "Staring")
        display(p)
        fileName = "/home/wesseling/DataDisk/Wesseling/WesW/Projects/Achilles/Plots/kritiek_" * string(aBofekId) * "_" * string(Int(10.0*aFlux)) * "_" * string(aNumber) * ".png"
	      savefig(p, fileName)
      catch ex
        println("???ERROR in plotCriticalValues: ", ex)
      end
    finally
    end
  end

  function calculate()
    global bofekProfile = DataManager.readBofekProfiles()
    global staring = DataManager.readStaringParams()
    global flux = DataManager.readFluxes()

    gr()

    what = "Ksat"
    if what == "KritiekeZ"
      q = [0.1, 0.2]
      n = 10000
      if test
      	n = 100
      end
      for i in 1:size(bofekProfile,1)
        if bofekProfile[i].bofekId == 2003
     #     bofekProfile[i].bofekId == 3005 ||
     #     bofekProfile[i].bofekId == 3017 ||
     #     bofekProfile[i].bofekId == 4003 ||
     #     bofekProfile[i].bofekId == 4007 ||
     #     bofekProfile[i].bofekId == 4009 ||
     #     bofekProfile[i].bofekId == 4014 ||
     #     bofekProfile[i].bofekId == 4017
          for j in 1:size(q,1)
            println(i, " ",j)
#           Old samples
            period = "old"
            global soilLayer = DataManager.readSoilLayers(bofekProfile[i].bofekId)
            createArrayWithParams(soilLayer, period)
            println("Number of combinations:",size(combination,1))
            generateZValues(q[j], n, period)
#           New samples
            period = "new"
            global soilLayer = DataManager.readSoilLayers(bofekProfile[i].bofekId)
            createArrayWithParams(soilLayer, period)
            println("Number of combinations:",size(combination,1))
            if size(combination,1) > 1
              generateZValues(q[j], n, period)
            end
#           Staring
            for i in 1:size(soilLayer,1)
              soilLayer[i].param = staring[soilLayer[i].staringId]
            end
            myValues = computeCriticalValues(q[j])
            println(myValues)
            plotCriticalValues(myValues, bofekProfile[i].bofekId, q[j], n)
          end
        end
      end
    end

    if what == "Ksat"
      exploreKsat()
    end
    if what == "Swap"
      location = "kat"  # kat or nba
      testWithSwap(location)
    end
    if what == "Storage"
      period = "new"
      n = 10000
      for i in 1:size(bofekProfile,1)
        if bofekProfile[i].bofekId == 2003 ||
            bofekProfile[i].bofekId == 3005 ||
            bofekProfile[i].bofekId == 3017 ||
            bofekProfile[i].bofekId == 4003 ||
            bofekProfile[i].bofekId == 4007 ||
            bofekProfile[i].bofekId == 4009 ||
            bofekProfile[i].bofekId == 4014 ||
            bofekProfile[i].bofekId == 4017
          for j in 1:size(flux,1)
            println(i, " ",j)
            global soilLayer = DataManager.readSoilLayers(bofekProfile[i].bofekId)
#          println(size(soilLayer,1), " layers")
            createArrayWithParams(soilLayer, period)
            println("Number of combinations:",size(combination,1))
            if size(combination,1) > 0
              generateValues(flux[j], n, period)
              computeWithStaring(flux[j])
              plotStorage(bofekProfile[i].bofekId, flux[j], n, period)
            end
          end
        end
      end
    end
  end

end

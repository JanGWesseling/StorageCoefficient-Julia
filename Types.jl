module Types

  mutable struct Flux
    id :: Int64
    value :: Float64
    name :: String
  end

  mutable struct SwapOutput
    gwlMin :: Float64
    gwlMax :: Float64
    gwlAve :: Float64
    qDrain :: Float64
    qDeep :: Float64
    epp :: Float64
    epa :: Float64
    esp :: Float64
    esa :: Float64
    potato :: Float64
    grass :: Float64
  end

  mutable struct SwapOutput2
    gwlMin :: Float64
    gwlMax :: Float64
    gwlAve :: Float64
    qDrain :: Float64
    qDeep :: Float64
    epp :: Float64
    epa :: Float64
    esp :: Float64
    esa :: Float64
    potato :: Float64
    pMin :: Float64
    pMax :: Float64
    t10 :: Int64
    t20 :: Int64
    t30 :: Int64
    t40 :: Int64
    t50 :: Int64
    t75 :: Int64
    t100 :: Int64
  end

  mutable struct BofekProfile
    bofekId :: Int64
    profileId :: Int64
    soilType :: String
  end

  mutable struct MVGParam
    id :: Int64
    thetaR :: Float64
    thetaS :: Float64
    kSat :: Float64
    alpha :: Float64
    l :: Float64
    n :: Float64
    m :: Float64
  end

  mutable struct SoilLayer
    id :: Int64
    size :: Int64
    top :: Int64
    bottom :: Int64
    staringId :: Int64
    param :: MVGParam
  end

  mutable struct Sample
    param :: Array{Main.Control.Types.MVGParam}
  end

  mutable struct Critical
    distance :: Float64
    volume :: Float64
  end

end

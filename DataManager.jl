module DataManager

  using MySQL
  using DataFrames

  function readFluxes()
    flux = Array{Main.Control.Types.Flux}(undef,0)
    try
      try
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sql = """SELECT * FROM fluxes;"""
        result = DBInterface.execute(myConnection, sql) |> DataFrame
        n = size(result,1)
        resize!(flux, n)
        for i in 1:n
          myFlux = Main.Control.Types.Flux(result[i,1], result[i,2], result[i,3])
          flux[i] = myFlux
        end
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in readFluxes: ",ex)
      end
    finally
    end
    return flux
  end

  function readBofekProfiles()
    profile = Array{Main.Control.Types.BofekProfile}(undef,0)
    try
      try
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sql = """SELECT * FROM bofekprofiles;"""
        result = DBInterface.execute(myConnection, sql) |> DataFrame
        n = size(result,1)
        resize!(profile, n)
        for i in 1:n
          myProfile = Main.Control.Types.BofekProfile(result[i,1], result[i,2], result[i,3])
          profile[i] = myProfile
        end
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in readBofekProfiles: ",ex)
      end
    finally
    end
    return profile
  end

  function readSoilLayers(aBofekId :: Int64)
    layer = Array{Main.Control.Types.SoilLayer}(undef,0)
    try
      try
        myParams = Main.Control.Types.MVGParam(-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sql = "SELECT * FROM soilprofiles WHERE bofekid=" * string(aBofekId) * " ORDER BY layerid;"
        result = DBInterface.execute(myConnection, sql) |> DataFrame
        n = size(result,1)
        resize!(layer, n)
        top = 0.0
        bottom = 0.0
        for i in 1:n
          top = bottom
          bottom = bottom - result[i,3]
#         extend profile to 5 m
          if i==n
            bottom = -310.0
          end
          myLayer = Main.Control.Types.SoilLayer(result[i,2], result[i,3], top, bottom, result[i,4], myParams)
          layer[i] = myLayer
        end
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in readSoilLayers: ",ex)
      end
    finally
    end
    return layer
  end

  function readStaringParams()
    staring = Array{Main.Control.Types.MVGParam}(undef,0)
    try
      try
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sql = "SELECT * FROM staringparams ORDER BY staringid;"
        result = DBInterface.execute(myConnection, sql) |> DataFrame
        n = size(result,1)
        resize!(staring, n)
        for i in 1:n
          m = 1.0 - (1.0 / result[i,7])
          myStaring = Main.Control.Types.MVGParam(result[i,1], result[i,2], result[i,3], result[i,4], result[i,5], result[i,6], result[i,7], m)
          staring[i] = myStaring
        end
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in readStaringParams: ",ex)
      end
    finally
    end
    return staring
  end

  function readSampleData(aStaringId :: Int64, aPeriod :: String)
    samples = Array{Main.Control.Types.MVGParam}(undef,0)
    try
      try
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sql = ""
        if aPeriod == "all"
          sql = "SELECT * FROM samples WHERE staringid=" * string(aStaringId) * " ORDER BY number;"
        end
        if aPeriod == "old"
          sql = "SELECT * FROM samples WHERE staringid=" * string(aStaringId) * " AND number < 10000 ORDER BY number;"
        end
        if aPeriod == "new"
          sql = "SELECT * FROM samples WHERE staringid=" * string(aStaringId) * " AND number >= 10000 ORDER BY number;"
        end

        result = DBInterface.execute(myConnection, sql) |> DataFrame
        n = size(result,1)
  #      println("nr of samples = ",n)
        resize!(samples, n)
        for i in 1:n
          mySample = Main.Control.Types.MVGParam(result[i,2], result[i,3], result[i,4], result[i,5], result[i,6], result[i,7], result[i,8], result[i,9])
          samples[i] = mySample
        end
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in readSampleData: ",ex)
      end
    finally
    end
    return samples
  end

  function readSampleData2(aId :: Int64)
    samples = Array{Main.Control.Types.MVGParam}(undef,0)
    try
      try
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sql = "SELECT * FROM samples WHERE number=" * string(aId) * ";"
        result = DBInterface.execute(myConnection, sql) |> DataFrame
        n = size(result,1)
        resize!(samples, n)
        for i in 1:n
          mySample = Main.Control.Types.MVGParam(result[i,2], result[i,3], result[i,4], result[i,5], result[i,6], result[i,7], result[i,8], result[i,9])
          samples[i] = mySample
        end
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in readSampleData2: ",ex)
      end
    finally
    end
    return samples
  end

  function deleteFromNBA()
    try
      try
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sql = "DELETE FROM nba;"
        result = DBInterface.execute(myConnection, sql)
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in deleteFromNBA: ",ex)
      end
    finally
    end
  end

  function storeSwapData(aRunId :: Int64, aData :: Main.Control.Types.SwapOutput2, aLocation :: String)
    try
      try
        myConnection = DBInterface.connect(MySQL.Connection,"127.0.0.1", "bofek","bofek", db="bofek")
        sqlDelete = "DELETE FROM " * aLocation * " WHERE runid=" * string(aRunId) * ";"
        result = DBInterface.execute(myConnection, sqlDelete)
        sqlInsert = "INSERT INTO " * aLocation * " (runid, gwlmin, gwlmax, gwlave, qdrain, qdeep, epp, epa, esp, esa, potato, pmin, pmax, t10, t20, t30, t40, t50, t75, t100) VALUE (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
        myStatement = DBInterface.prepare(myConnection, sqlInsert)
        result = DBInterface.execute(myStatement,[aRunId, aData.gwlMin, aData.gwlMax, aData.gwlAve,
           aData.qDrain, aData.qDeep, aData.epp, aData.epa, aData.esp, aData.esa,
           aData.potato, aData.pMin, aData.pMax, aData.t10, aData.t20, aData.t30, aData.t40, aData.t50, aData.t75, aData.t100])
        DBInterface.close!(myConnection)
      catch ex
        println("???ERROR in storeSwapData: ",ex)
      end
    finally
    end
  end


end

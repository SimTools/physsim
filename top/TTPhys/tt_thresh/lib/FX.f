C
C  *INCLUDE F1
C  *INCLUDE F2
C        IMPLICIT REAL*8 ( A-H, O-Z )
C  C--
C        B0 = (33 -2*5)/6.D0
C        B1 = (153 -19*5)/12.D0
C        A    = B1/B0**2
C  C--
C        NX   = 100
C        XMIN =  11
C        XMAX = 1501
C        DX   = (XMAX-XMIN)/NX
C  C--
C        WRITE(20,*) 'NEW FRAME'
C        WRITE(20,*) 'SET LIMIT X ', XMIN, ' ', XMAX
C        WRITE(20,*) 'SET LIMIT Y 0 3'
C        WRITE(20,*) 'SET ORDER X Y'
C        DO 10 IX = 0, NX
C           X  = XMIN + IX*DX
C           Y  = F1(X) + A*F2(X)
C           WRITE(20,*) X, Y
C  10    CONTINUE
C        WRITE(20,*) 'JOIN DOT'
C  C--
C        NX   =  300
C        DX   = (XMAX-XMIN)/NX
C        DO 20 IX = 0, NX
C           X  = XMIN + IX*DX
C           Y  = FX(X)
C           WRITE(20,*) X, Y
C  20    CONTINUE
C        WRITE(20,*) 'JOIN'
C        END
C*
C*
C*
 
      DOUBLE PRECISION FUNCTION FX(X)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   X
C--
      REAL*8 FXDATA(0:300)
      DATA NX     /       300 /
      DATA XMIN   /.10000D+02 /
      DATA DX     /.50000D+01 /
      DATA (FXDATA(K),K=  0, 59) /
     . 0.1331D+01, 0.1385D+01, 0.1442D+01, 0.1540D+01, 0.1570D+01,
     . 0.1582D+01, 0.1631D+01, 0.1662D+01, 0.1664D+01, 0.1688D+01,
     . 0.1719D+01, 0.1723D+01, 0.1732D+01, 0.1757D+01, 0.1767D+01,
     . 0.1768D+01, 0.1785D+01, 0.1801D+01, 0.1801D+01, 0.1809D+01,
     . 0.1825D+01, 0.1829D+01, 0.1831D+01, 0.1844D+01, 0.1852D+01,
     . 0.1852D+01, 0.1860D+01, 0.1871D+01, 0.1872D+01, 0.1874D+01,
     . 0.1885D+01, 0.1890D+01, 0.1890D+01, 0.1897D+01, 0.1905D+01,
     . 0.1905D+01, 0.1908D+01, 0.1917D+01, 0.1919D+01, 0.1919D+01,
     . 0.1926D+01, 0.1932D+01, 0.1931D+01, 0.1935D+01, 0.1942D+01,
     . 0.1943D+01, 0.1944D+01, 0.1951D+01, 0.1954D+01, 0.1953D+01,
     . 0.1958D+01, 0.1963D+01, 0.1963D+01, 0.1965D+01, 0.1971D+01,
     . 0.1973D+01, 0.1972D+01, 0.1977D+01, 0.1981D+01, 0.1981D+01/
      DATA (FXDATA(K),K= 60,119) /
     . 0.1983D+01, 0.1988D+01, 0.1989D+01, 0.1989D+01, 0.1994D+01,
     . 0.1997D+01, 0.1996D+01, 0.1999D+01, 0.2003D+01, 0.2003D+01,
     . 0.2004D+01, 0.2008D+01, 0.2010D+01, 0.2010D+01, 0.2013D+01,
     . 0.2016D+01, 0.2016D+01, 0.2018D+01, 0.2022D+01, 0.2022D+01,
     . 0.2022D+01, 0.2026D+01, 0.2028D+01, 0.2028D+01, 0.2030D+01,
     . 0.2033D+01, 0.2033D+01, 0.2034D+01, 0.2037D+01, 0.2039D+01,
     . 0.2039D+01, 0.2041D+01, 0.2044D+01, 0.2044D+01, 0.2044D+01,
     . 0.2048D+01, 0.2049D+01, 0.2048D+01, 0.2051D+01, 0.2053D+01,
     . 0.2053D+01, 0.2054D+01, 0.2057D+01, 0.2058D+01, 0.2058D+01,
     . 0.2061D+01, 0.2062D+01, 0.2062D+01, 0.2063D+01, 0.2066D+01,
     . 0.2066D+01, 0.2066D+01, 0.2069D+01, 0.2070D+01, 0.2070D+01,
     . 0.2072D+01, 0.2074D+01, 0.2074D+01, 0.2075D+01, 0.2077D+01/
      DATA (FXDATA(K),K=120,179) /
     . 0.2078D+01, 0.2078D+01, 0.2080D+01, 0.2081D+01, 0.2081D+01,
     . 0.2082D+01, 0.2085D+01, 0.2085D+01, 0.2085D+01, 0.2087D+01,
     . 0.2088D+01, 0.2088D+01, 0.2090D+01, 0.2091D+01, 0.2091D+01,
     . 0.2092D+01, 0.2094D+01, 0.2095D+01, 0.2095D+01, 0.2096D+01,
     . 0.2098D+01, 0.2098D+01, 0.2099D+01, 0.2101D+01, 0.2101D+01,
     . 0.2101D+01, 0.2103D+01, 0.2104D+01, 0.2104D+01, 0.2105D+01,
     . 0.2107D+01, 0.2107D+01, 0.2107D+01, 0.2109D+01, 0.2110D+01,
     . 0.2109D+01, 0.2111D+01, 0.2112D+01, 0.2112D+01, 0.2113D+01,
     . 0.2115D+01, 0.2115D+01, 0.2115D+01, 0.2117D+01, 0.2118D+01,
     . 0.2117D+01, 0.2118D+01, 0.2120D+01, 0.2120D+01, 0.2120D+01,
     . 0.2122D+01, 0.2123D+01, 0.2122D+01, 0.2124D+01, 0.2125D+01,
     . 0.2125D+01, 0.2125D+01, 0.2127D+01, 0.2127D+01, 0.2127D+01/
      DATA (FXDATA(K),K=180,239) /
     . 0.2129D+01, 0.2130D+01, 0.2130D+01, 0.2130D+01, 0.2132D+01,
     . 0.2132D+01, 0.2132D+01, 0.2134D+01, 0.2134D+01, 0.2134D+01,
     . 0.2135D+01, 0.2136D+01, 0.2136D+01, 0.2137D+01, 0.2138D+01,
     . 0.2139D+01, 0.2138D+01, 0.2140D+01, 0.2141D+01, 0.2141D+01,
     . 0.2141D+01, 0.2143D+01, 0.2143D+01, 0.2143D+01, 0.2144D+01,
     . 0.2145D+01, 0.2145D+01, 0.2145D+01, 0.2147D+01, 0.2147D+01,
     . 0.2147D+01, 0.2148D+01, 0.2149D+01, 0.2149D+01, 0.2150D+01,
     . 0.2151D+01, 0.2151D+01, 0.2151D+01, 0.2152D+01, 0.2153D+01,
     . 0.2152D+01, 0.2154D+01, 0.2155D+01, 0.2154D+01, 0.2155D+01,
     . 0.2156D+01, 0.2156D+01, 0.2156D+01, 0.2158D+01, 0.2158D+01,
     . 0.2158D+01, 0.2159D+01, 0.2160D+01, 0.2160D+01, 0.2160D+01,
     . 0.2161D+01, 0.2162D+01, 0.2161D+01, 0.2162D+01, 0.2163D+01/
      DATA (FXDATA(K),K=240,300) /
     . 0.2163D+01, 0.2164D+01, 0.2165D+01, 0.2165D+01, 0.2165D+01,
     . 0.2166D+01, 0.2167D+01, 0.2166D+01, 0.2167D+01, 0.2168D+01,
     . 0.2168D+01, 0.2168D+01, 0.2169D+01, 0.2170D+01, 0.2170D+01,
     . 0.2171D+01, 0.2171D+01, 0.2171D+01, 0.2172D+01, 0.2173D+01,
     . 0.2173D+01, 0.2173D+01, 0.2174D+01, 0.2175D+01, 0.2174D+01,
     . 0.2175D+01, 0.2176D+01, 0.2176D+01, 0.2176D+01, 0.2177D+01,
     . 0.2178D+01, 0.2177D+01, 0.2178D+01, 0.2179D+01, 0.2179D+01,
     . 0.2179D+01, 0.2180D+01, 0.2180D+01, 0.2180D+01, 0.2181D+01,
     . 0.2182D+01, 0.2182D+01, 0.2182D+01, 0.2183D+01, 0.2183D+01,
     . 0.2183D+01, 0.2184D+01, 0.2185D+01, 0.2185D+01, 0.2185D+01,
     . 0.2186D+01, 0.2186D+01, 0.2186D+01, 0.2187D+01, 0.2187D+01,
     . 0.2187D+01, 0.2188D+01, 0.2189D+01, 0.2189D+01, 0.2189D+01,
     . 0.2190D+01/
C
C========< Entry Point >================================================
C
C--
C  Calculate f_x(x) = f_1(x) + b_1/b_0^2*f_2(x).
C--
      IX = (X-XMIN)/DX
      IF ( IX.LT.0 .OR. IX.GE.NX ) THEN
         PRINT *, ' Error in FX: X = ', X, ' out of range.'
      ELSE
         FX = FXDATA(IX) + (FXDATA(IX+1)-FXDATA(IX))*(X-XMIN-IX*DX)/DX
      END IF
C--
C  That's it.
C--
      RETURN
      END

module thermalConstants

    implicit none

    real, parameter :: airHeatCap=1.006d3, tempAir=25.d0+273.d0, tempAir4=tempAir**4, lw=2256.d3
    real            :: skinDensityInit

    contains

        real function airThermalCond(T, i)
        !data taken from engineeringtoolbox.com and fitted in gnuplot
        !temp input in kelvin
            implicit none

            real, intent(IN) :: T!in Kelvin
            real, parameter  :: a=-0.188521, b=0.000367259, c=0.212453
            integer :: i

            if(t < 0.)then
                print*,i,'loop #'
                error stop "Negative temp in kelvin"
            end if
            airThermalCond = a * exp(-b * (T - 273.15d0)) + c
        end function airThermalCond


        real function airDensity(T)
        !values taken from wikipedia
        !temp input in kelvin
            implicit none

            real, intent(IN) :: T !in kelvin
            real, parameter  :: pressue1ATM=101.325d3, Rspec=287.058d0

            airDensity = pressue1ATM / (Rspec * T)

        end function airDensity


        !source for three below function: Analysis of Thermal Relaxation During LaserIrradiation of Tissue, Choi & Welch 2001
        real function getSkinDensity()
        !get dsensity of skin based upon water content. In Kg m-3

            implicit none

            getSkinDensity = 1120.d0

        end function getSkinDensity


        real function getSkinHeatCap()
        !get heat capacity for skin. In J Kg-1 K-1 

            implicit none

            getSkinHeatCap = 3200.d0

        end function getSkinHeatCap


        real function getSkinThermalCond()
        !get thermal conductivity for skin. In W m-1 K-1 

            implicit none

            getSkinThermalCond = .34d0!currentDensity * (6.28d-4*waterContent + 1.17d-4*proteinContent)!new function with protein

        end function getSkinThermalCond
end module thermalConstants
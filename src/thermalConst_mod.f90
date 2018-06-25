module thermalConstants

    implicit none

    real, parameter :: skinThermalCond=0.209d0, skinDensity=1070.d0, skinHeatCap=3400.d0
    real, parameter :: airHeatCap=1.006d3, tempAir=25.d0+273.d0, tempAir4=tempAir**4
    real, parameter :: S_B_Constant=5.670367d-8, h=10.d0, eps=0.98d0,lw=2256.d3
    real, parameter :: skinAlpha=skinThermalCond/(skinDensity * skinHeatCap)
    real            :: skinBeta,QVapor

    !skinDensity    : 1020 Kg/m3    source 10.1109/TBME.2009.2033464  burn papers says: 1093-1190
    !skinThermalCond: 0.21 W/m K    source burn paper- pig skin
    !skinHeatCap    : 3349 J/Kg K   source burn paper- pig skin

    !airHeatCap: J/kg K  source engineeringtoolbox.com 

    !lw: J/kg latent heat of vaporisation water

    ! S_B_Constant: stefan-Boltzmann constant
    ! h: convective constant W/m2 K  between 2-10                   ************need source************
    ! eps: emissivity                           ************need source************

    contains

        real function airThermalCond(T)
        !data taken from engineeringtoolbox.com and fitted in gnuplot
            implicit none

            real, intent(IN) :: T!in Kelvin
            real, parameter  :: a=-0.188521, b=0.000367259, c=0.212453

            airThermalCond = a * exp(-b * (T - 273.15d0)) + c

        end function airThermalCond


        real function airDensity(T)
        !values taken from wikipedia
            implicit none

            real, intent(IN) :: T !in kelvin
            real, parameter  :: pressue1ATM=101.325d3, Rspec=287.058d0

            airDensity = pressue1ATM / (Rspec * T)

        end function airDensity
end module thermalConstants

        !init heat variables for medium
        ! kappa = 0.00209 !0.0056 ! W/cm K    thermal conductivity
        ! rho = 1.07 ! g/cm^3                 density
        ! c_heat = 3.4 !J/g K                 heat capacity
        ! eps = 0.98
        ! sigma = 5.670373e-8
        ! t_air = 25.+273.
        ! t_air4 = t_air**4
        ! h = 10.
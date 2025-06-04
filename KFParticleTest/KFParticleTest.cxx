#include <cmath>

#include "KFParticle.hxx"
#include "KFParticle_Math.hxx"

int main() {
    constexpr double pion_mass{0.1395704};  // GeV/c
    constexpr double k0s_mass{0.49761100};  // GeV/c
    constexpr double b_field{4.9798};

    // Vertex 0 //

    KF::Vector<3> vv0{0., 0., 10.};
    KF::SymMatrix<3> cc0{0.01, 0., 0.01, 0., 0., 0.01};

    // Particle 1 //

    // clang-format off
    KF::Vector<6> pp1{-0.061996019110347252, -1.3579236865955473,   27.147283554077148,
                      0.62539337626870062,   -0.028552340672283318, -0.18467358509984011};
    KF::SymMatrix<6> cc1{3.3055800809774214E-05,
                         0.00098316976438185002,  0.04740889543423539,
                         -8.5596097466772512E-05, -0.0037516094381694971,  0.032156504690647125,
                         -2.2812597903705375E-05, -0.0012121012247057524,  3.0646383360925928E-05,  6.1388628418184652E-05,
                         -4.4071909055788304E-06, -0.00048870318030618627, 3.8062554692505919E-05,  1.2177141510445709E-05, 7.6900178535210476E-06,
                         6.6224441962932268E-06,  0.00034363110217286891,  -0.00031520420397528146, -1.6277704753223909E-05, -3.4322154557097545E-06, 1.027411488502718E-05};
    // clang-format on
    // track1.SetNDF(1);
    // track1.SetChi2(1.5);
    KF::Particle p1{pp1, cc1, -1, pion_mass};

    // Particle 2 //

    // clang-format off
    KF::Vector<6> pp2{-0.20371287092090862, 3.0678058943547839,   -19.93988037109375,
                      0.37533048135363339,  0.024923235867488316, 0.19031024520542122};
    KF::SymMatrix<6> cc2{0.00022312908970259721,
                         -0.00064291160449645151, 0.089331037457232143,
                         0.00047880877483649206,  -0.045478494677353445,   0.11199165135622025,
                         4.6362085390124077E-07,  0.00070978326424729935,  -0.00014164977426380486, 1.7553871209443515E-05,
                         -2.2044831998838091E-05, -0.00059994741249631909, 0.00030148707952079015,  -4.6574515272730461E-06, 7.2618497455845866E-06,
                         -1.2427988441207971E-06, 0.00030830063771211896,  -0.00061853865528922161, 5.4390968700069889E-06, -1.9914477627292868E-06, 8.9837108094398403E-06};
    // clang-format on
    // track2.SetNDF(2);
    // track2.SetChi2(2.5);
    KF::Particle p2{pp2, cc2, +1, pion_mass};

    // Composite Particle A (mass constraint OFF, vertex hypothesis OFF) //

    std::cout << "===================================================================================" << '\n';
    std::cout << "== START Composite Particle A (mass constraint OFF, vertex hypothesis OFF) START ==" << '\n';
    std::cout << "===================================================================================" << '\n';
    KF::Particle cpA;
    cpA.AddDaughter(p1, b_field);
    cpA.AddDaughter(p2, b_field);
    cpA.Print();
    std::cout << "=================================================================================" << '\n';
    std::cout << "== END Composite Particle A (mass constraint OFF, vertex hypothesis OFF) START ==" << '\n';
    std::cout << "=================================================================================" << '\n';

    std::cout << '\n';

    // Composite Particle B (mass constraint ON, vertex hypothesis OFF) //

    std::cout << "==================================================================================" << '\n';
    std::cout << "== START Composite Particle B (mass constraint ON, vertex hypothesis OFF) START ==" << '\n';
    std::cout << "==================================================================================" << '\n';
    KF::Particle cpB;
    cpB.AddDaughter(p1, b_field);
    cpB.AddDaughter(p2, b_field);
    cpB.AddMassConstraint(k0s_mass);
    cpB.Print();
    std::cout << "==============================================================================" << '\n';
    std::cout << "== END Composite Particle B (mass constraint ON, vertex hypothesis OFF) END ==" << '\n';
    std::cout << "==============================================================================" << '\n';

    std::cout << '\n';

    // Composite Particle C (mass constraint OFF, vertex hypothesis ON) //

    std::cout << "==================================================================================" << '\n';
    std::cout << "== START Composite Particle C (mass constraint OFF, vertex hypothesis ON) START ==" << '\n';
    std::cout << "==================================================================================" << '\n';
    KF::Particle cpC;
    cpC.AddDaughter(p1, b_field);
    cpC.AddDaughter(p2, b_field);
    cpC.AddProductionVertex(vv0, cc0, b_field);
    cpC.Print();
    std::cout << "==============================================================================" << '\n';
    std::cout << "== END Composite Particle C (mass constraint OFF, vertex hypothesis ON) END ==" << '\n';
    std::cout << "==============================================================================" << '\n';

    std::cout << '\n';

    // Composite Particle D (mass constraint ON, vertex hypothesis ON) //

    std::cout << "=================================================================================" << '\n';
    std::cout << "== START Composite Particle D (mass constraint ON, vertex hypothesis ON) START ==" << '\n';
    std::cout << "=================================================================================" << '\n';
    KF::Particle cpD;
    cpD.AddDaughter(p1, b_field);
    cpD.AddDaughter(p2, b_field);
    cpD.AddMassConstraint(k0s_mass);
    cpD.AddProductionVertex(vv0, cc0, b_field);
    cpD.Print();
    std::cout << "=============================================================================" << '\n';
    std::cout << "== END Composite Particle D (mass constraint ON, vertex hypothesis ON) END ==" << '\n';
    std::cout << "=============================================================================" << '\n';
}

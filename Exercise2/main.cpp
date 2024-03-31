#include <iostream>
#include "Eigen/Eigen"

using namespace Eigen;
using namespace std;

// Funzione per calcolare Errore
double ErroreRelativo(const VectorXd& xSoluzione, const VectorXd& xAtteso) {
    return (xSoluzione - xAtteso).norm() / xAtteso.norm();
}

int main() {

    Matrix2d A1, A2, A3;
    Vector2d b1, b2, b3;

    A1 << 5.547001962252291e-01, -3.770900990025203e-02,
        8.320502943378437e-01, -9.992887623566787e-01;
    b1 << -5.169911863249772e-01, 1.672384680188350e-01;

    A2 << 5.547001962252291e-01, -5.540607316466765e-01,
        8.320502943378437e-01, -8.324762492991313e-01;
    b2 << -6.394645785530173e-04, 4.259549612877223e-04;

    A3 << 5.547001962252291e-01, -5.547001955851905e-01,
        8.320502943378437e-01, -8.320502947645361e-01;
    b3 << -6.400391328043042e-10, 4.266924591433963e-10;

    // soluzione dei sistemi utilizzando la fattorizzazione PALU
    Vector2d x1_PALU = A1.partialPivLu().solve(b1);
    Vector2d x2_PALU = A2.partialPivLu().solve(b2);
    Vector2d x3_PALU = A3.partialPivLu().solve(b3);

    // soluzione utilizzando la fattorizzazione QR
    Vector2d x1_QR = A1.colPivHouseholderQr().solve(b1);
    Vector2d x2_QR = A2.colPivHouseholderQr().solve(b2);
    Vector2d x3_QR = A3.colPivHouseholderQr().solve(b3);


    Vector2d x_atteso(-1.0e+0, -1.0e+0);

    // Calcolo errori relativi i sistemi
    double error1_PALU = ErroreRelativo(x1_PALU, x_atteso);
    double error2_PALU = ErroreRelativo(x2_PALU, x_atteso);
    double error3_PALU = ErroreRelativo(x3_PALU, x_atteso);

    double error1_QR = ErroreRelativo(x1_QR, x_atteso);
    double error2_QR = ErroreRelativo(x2_QR, x_atteso);
    double error3_QR = ErroreRelativo(x3_QR, x_atteso);


    cout << "\nSistema 1" << endl;
    cout << "Soluzione PALU: " << x1_PALU.transpose() << ", Errore relativo: " << error1_PALU << endl;
    cout << "Soluzione QR: " << x1_QR.transpose() << ", Errore relativo: " << error1_QR << endl;

    cout << "\nSistema 2:" << endl;
    cout << "Soluzione PALU: " << x2_PALU.transpose() << ", Errore relativo: " << error2_PALU << endl;
    cout << "Soluzione QR: " << x2_QR.transpose() << ", Errore relativo: " << error2_QR << endl;

    cout << "\nSistema 3:" << endl;
    cout << "Soluzione PALU: " << x3_PALU.transpose() << ", Errore relativo: " << error3_PALU << endl;
    cout << "Soluzione QR: " << x3_QR.transpose() << ", Errore relativo: " << error3_QR << endl;

    return 0;
}

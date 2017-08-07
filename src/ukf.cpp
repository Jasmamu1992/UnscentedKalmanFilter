#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // weights_computed is set to false
  weights_computed = false;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_aug){
    n_x_ = x_.size();
    n_aug_ = n_x_ + 2;
    lambda_ = 3 - n_aug_;

    //create augmented mean state
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_, n_x_) = std_a_ * std_a_;
    P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    //MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug->col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++){
        Xsig_aug->col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug->col(n_aug_ + i + 1) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }
}

double CircularPhi(double Phi){
    while(Phi < -M_PI) Phi += 2 * M_PI;
    while(Phi >  M_PI) Phi -= 2 * M_PI;
    return Phi;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    n_x_ = x_.size();
    n_aug_ = n_x_ + 2;
    MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    AugmentedSigmaPoints(&Xsig_aug_);
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        VectorXd xaug = Xsig_aug_.col(i);
        double p_x = xaug(0);
        double p_y = xaug(1);
        double v   = xaug(2);
        double yaw = xaug(3);
        double yawd = xaug(4);
        double nu_a = xaug(5);
        double nu_yawdd = xaug(6);
        double px_p;
        double py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //adding noise
        px_p += 0.5 * delta_t * delta_t * cos(yaw) * nu_a;
        py_p += 0.5 * delta_t * delta_t * sin(yaw) * nu_a;
        v_p += delta_t * nu_a;
        yaw_p += 0.5 * delta_t * delta_t * nu_yawdd;
        yawd_p += delta_t * nu_yawdd;

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if (!is_initialized_){
        is_initialized_ = true;
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
            float ro = meas_package.raw_measurements_(0);
            float theta = meas_package.raw_measurements_(1);
            float ro_dot = meas_package.raw_measurements_(2);
            float Px = ro * cos(theta);
            float Py = ro * sin(theta);
            float V = ro_dot;
            x_ << Px, Py, V, 0, 0;
            P_ = MatrixXd::Identity(x_.size(), x_.size());
            P_ *= 0.5;
            time_us_ = meas_package.timestamp_;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
            float Px = meas_package.raw_measurements_(0);
            float Py = meas_package.raw_measurements_(1);
            x_ << Px, Py, 0, 0, 0;
            P_ = MatrixXd::Identity(x_.size(), x_.size());
            P_ *= 0.5;
            time_us_ = meas_package.timestamp_;
        }
        return;
    }
    else{
        float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
        time_us_ = meas_package.timestamp_;
        Prediction(delta_t);
        if (meas_package.sensor_type_ == MeasurementPackage::LASER){
            UpdateLidar(meas_package);
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
            UpdateRadar(meas_package);
        }
    }
}

void UKF::ComputeWeights(){
    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    float weight = 1 / (2 * (lambda_ + n_aug_));
    for (int i = 1; i < 2 * n_aug_ + 1; i++){
        weights_(i) = weight;
    }
}

void UKF::PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred, MatrixXd Xsig, int AngleInd){
    if (!weights_computed){
        ComputeWeights();
        weights_computed = true;
    }
    int n = x_pred->size();
    VectorXd x = VectorXd(n);
    MatrixXd P = MatrixXd(n, n);
    x.fill(0);
    P.fill(0);

    //compute mean
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        x = x + Xsig.col(i) * weights_(i);
    }

    //compute covariance
    for (int i = 1; i < 2 * n_aug_ + 1; i++){
        VectorXd x_diff = Xsig.col(i) - x;

        x_diff(AngleInd) = CircularPhi(x_diff(AngleInd));

        P = P + weights_(i) * (Xsig.col(i) - x) * (Xsig.col(i) - x).transpose();
    }

    *x_pred = x;
    *P_pred = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    MatrixXd P_pred = MatrixXd(n_x_, n_x_);
    VectorXd x_pred = VectorXd(n_x_);

    PredictMeanAndCovariance(&x_pred, &P_pred, Xsig_pred_, 3);
    x_ = x_pred;
    P_ = P_pred;

    int nl = meas_package.raw_measurements_.size();
    VectorXd z = meas_package.raw_measurements_;

    MatrixXd H_(nl, n_x_);
    H_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;

    MatrixXd R_(nl, nl);
    R_ << std_laspx_*std_laspx_,     0,
          0,                         std_laspy_*std_laspy_;

    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;

    x_(3) = CircularPhi(x_(3));
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
prediction
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and cov3ariance, P_.

  You'll also need to calculate the radar NIS.
  */
    VectorXd z = meas_package.raw_measurements_;
    int nr = z.size();
    MatrixXd R = MatrixXd(nr, nr);
    R << std_radr_ * std_radr_, 0,                         0,
         0,                     std_radphi_ * std_radphi_, 0,
         0,                     0,                         std_radrd_ * std_radrd_;
    MatrixXd Zsig(nr, 2 * n_aug_ + 1);

    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        //yaw = CircularPhi(yaw);
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        //Zsig(1,i) = CircularPhi(Zsig(1,i));
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);                 //rd
    }

    //Calculate Mean of predicted sigma points in polar domain
    MatrixXd P_pred = MatrixXd(nr, nr);
    VectorXd z_pred = VectorXd(nr);
    PredictMeanAndCovariance(&z_pred, &P_pred, Zsig, 1);

    //Calculate Mean of predicted sigma points in cartesian domain
    MatrixXd P = MatrixXd(n_x_, n_x_);
    VectorXd x_pred = VectorXd(n_x_);
    PredictMeanAndCovariance(&x_pred, &P, Xsig_pred_, 3);
    x_ = x_pred;
    P_ = P;

    //Calculating S
    MatrixXd S = P_pred + R;

    //Calculating Tc
    MatrixXd Tc = MatrixXd(n_x_, nr);
    Tc.fill(0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = CircularPhi(z_diff(1));

        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = CircularPhi(x_diff(3));

        if (!weights_computed){
            ComputeWeights();
            weights_computed = true;
        }
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;
    z_diff(1) = CircularPhi(z_diff(1));

    //Update State mean and Covariance Matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    x_(3) = CircularPhi(x_(3));

}

/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

#include "material_drucker_prager.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
/// Deviatoric Stress

/* -------------------------------------------------------------------------- */
/// Yield Stress
template <Int dim>
inline Real
MaterialDruckerPrager<dim>::computeYieldStress(const Matrix<Real> & sigma) {
  return this->alpha * sigma.trace() - this->k;
}

/* -------------------------------------------------------------------------- */
/// Yield function
template <Int dim>
inline Real
MaterialDruckerPrager<dim>::computeYieldFunction(const Matrix<Real> & sigma) {
  Matrix<Real, dim, dim> sigma_dev = Material::computeDeviatoric<dim>(sigma);

  // compute deviatoric invariant J2
  auto j2 = (1. / 2.) * sigma_dev.doubleDot(sigma_dev);
  auto sigma_dev_eff = std::sqrt(3. * j2);

  auto modified_yield_stress = computeYieldStress(sigma);

  return sigma_dev_eff + modified_yield_stress;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1, typename D2, typename D3,
          aka::enable_if_t<aka::are_vectors_v<D2, D3>> *>
inline void MaterialDruckerPrager<dim>::computeGradientAndPlasticMultplier(
    const Eigen::MatrixBase<D1> & sigma_trial, Real & plastic_multiplier_guess,
    Eigen::MatrixBase<D2> & gradient_f,
    Eigen::MatrixBase<D3> & delta_inelastic_strain, Int max_iterations,
    Real tolerance) {

  const Int size = voigt_h::size;

  // guess stress state at each iteration, initial guess is the trial state
  Matrix<Real, dim, dim> sigma_guess(sigma_trial);

  // plastic multiplier guess at each iteration, initial guess is zero
  plastic_multiplier_guess = 0.;

  // gradient of yield surface in voigt notation
  gradient_f.zero();

  // plastic strain increment at each iteration
  delta_inelastic_strain.zero();

  // variation in sigma at each iteration
  Vector<Real, size> delta_sigma;

  // krocker delta vector in voigt notation
  Vector<Real, size> kronecker_delta;
  kronecker_delta.template block<dim, 1>(0, 0).fill(1.);

  // hessian matrix of yield surface
  Matrix<Real, size, size> hessian_f;

  // scaling matrix for computing gradient and hessian from voigt notation
  Matrix<Real, size, size> scaling_matrix =
      Matrix<Real, size, size>::Identity();
  scaling_matrix.template block<dim, dim>(0, 0) =
      scaling_matrix.template block<dim, dim>(0, 0) * 2;

  // elastic stifnness tensor
  Matrix<Real, size, size> De;
  MaterialElastic<dim>::computeTangentModuliOnQuad(
      make_named_tuple("tangent_moduli"_n = De));

  // elastic compliance tensor
  Matrix<Real, size, size> Ce = De.inverse();

  // objective function to be computed
  Vector<Real, size> f;
  f.zero();

  // yield function value at each iteration
  Real yield_function;

  // if sigma is above the threshold value
  auto above_threshold = [&sigma_guess](Real & k, Real & alpha) {
    auto I1 = sigma_guess.trace();
    return I1 >= k / alpha;
  };

  // to project stress state at origin of yield function if first
  // invariant is greater than the threshold
  if (above_threshold(k, alpha) and this->alpha > 0) {
    auto update_first_obj = [&sigma_guess]() {
      auto sigma_dev = Material::computeDeviatoric<dim>(sigma_guess);
      auto error = (1. / 2) * sigma_dev.doubleDot(sigma_dev);
      return error;
    };

    auto update_sec_obj = [&sigma_guess](Real & k, Real & alpha) {
      auto error = alpha * sigma_guess.trace() - k;
      return error;
    };

    auto projection_error = update_first_obj();

    while (tolerance < projection_error) {
      auto sigma_dev = Material::computeDeviatoric<dim>(sigma_guess);

      sigma_guess += -projection_error * sigma_dev.inverse();
      projection_error = update_first_obj();
    }

    projection_error = update_sec_obj(k, alpha);

    while (tolerance < projection_error) {
      sigma_guess += -projection_error * 1. / this->alpha *
                     Matrix<Real, dim, dim>::Identity();
      projection_error = update_sec_obj(k, alpha);
    }

    auto delta_sigma_final = sigma_trial - sigma_guess;
    auto delta_sigma_voigt = voigt_h::matrixToVoigt(delta_sigma_final);

    delta_inelastic_strain = Ce * delta_sigma_voigt;

    return;
  }

  // lambda function to compute gradient of yield surface in voigt notation
  auto compute_gradient_f = [&sigma_guess, &scaling_matrix, &kronecker_delta,
                             &gradient_f](Real & alpha) {
    auto sigma_dev = Material::computeDeviatoric<dim>(sigma_guess);
    Vector<Real> sigma_dev_voigt = voigt_h::matrixToVoigt(sigma_dev);

    // compute deviatoric invariant
    auto j2 = (1. / 2.) * sigma_dev.doubleDot(sigma_dev);

    gradient_f =
        scaling_matrix * sigma_dev_voigt * 3. / (2. * std::sqrt(3. * j2)) +
        alpha * kronecker_delta;
  };

  // lambda function to compute hessian matrix of yield surface
  auto compute_hessian_f = [&sigma_guess, &hessian_f, &scaling_matrix,
                            &kronecker_delta]() {
    const Int size = voigt_h::size;
    auto sigma_dev = Material::computeDeviatoric<dim>(sigma_guess);
    auto sigma_dev_voigt = voigt_h::matrixToVoigt(sigma_dev);

    // compute deviatoric invariant J2
    auto j2 = (1. / 2.) * sigma_dev.doubleDot(sigma_dev);

    Vector<Real, size> temp = scaling_matrix * sigma_dev_voigt;

    auto id = -kronecker_delta * kronecker_delta.transpose() / 3. +
              Matrix<Real, size, size>::Identity();

    hessian_f =
        temp * temp.transpose() * -9. / (4. * std::pow(3. * j2, 3. / 2.)) +
        +(3. / (2. * std::pow(3. * j2, 1. / 2.))) * scaling_matrix * id;
  };

  /* --------------------------- */
  /* init before iteration loop  */
  /* --------------------------- */
  auto update_f = [&f, &sigma_guess, &sigma_trial, &plastic_multiplier_guess,
                   &Ce, &De, &yield_function, &gradient_f,
                   &delta_inelastic_strain,
                   &compute_gradient_f](Real & k, Real & alpha) {
    // compute gradient
    compute_gradient_f(alpha);

    // compute yield function
    auto sigma_dev = Material::computeDeviatoric<dim>(sigma_guess);
    auto j2 = (1. / 2.) * sigma_dev.doubleDot(sigma_dev);
    auto sigma_dev_eff = std::sqrt(3. * j2);

    auto modified_yield_stress = alpha * sigma_guess.trace() - k;

    yield_function = sigma_dev_eff + modified_yield_stress;

    // compute increment strain
    auto sigma_trial_voigt = voigt_h::matrixToVoigt(sigma_trial);
    auto sigma_guess_voigt = voigt_h::matrixToVoigt(sigma_guess);
    auto tmp = sigma_trial_voigt - sigma_guess_voigt;
    delta_inelastic_strain = Ce * tmp;

    // compute objective function
    f = De * gradient_f * plastic_multiplier_guess;
    f = tmp - f;

    // compute error
    auto error = std::max(f.norm(), std::abs(yield_function));
    return error;
  };

  Real alpha_tmp{alpha};
  Real k_tmp{k};
  if (radial_return_mapping) {
    alpha_tmp = 0;
    k_tmp = std::abs(alpha * sigma_guess.trace() - k);
  }

  auto projection_error = update_f(k_tmp, alpha_tmp);

  /* --------------------------- */
  /* iteration loop              */
  /* --------------------------- */
  Matrix<Real, size, size> xi;
  Matrix<Real, size, size> xi_inv;
  Vector<Real, size> tmp;
  Vector<Real, size> tmp1;
  Matrix<Real, size, size> tmp2;

  Int iterations{0};
  while (tolerance < projection_error and iterations < max_iterations) {

    // compute hessian at previous step
    compute_hessian_f();

    // compute inverse matrix Xi
    xi = Ce + plastic_multiplier_guess * hessian_f;

    xi_inv = xi.inverse();
    tmp = xi_inv * gradient_f;

    auto denominator = gradient_f.dot(tmp);

    // compute plastic multplier guess

    tmp1 = xi_inv * delta_inelastic_strain;
    plastic_multiplier_guess =
        (gradient_f.dot(tmp1) + yield_function) / denominator;

    // compute delta sigma

    tmp2 = xi_inv - tmp * tmp.transpose() / denominator;
    delta_sigma =
        tmp2 * delta_inelastic_strain - tmp * yield_function / denominator;

    // compute sigma_guess
    Matrix<Real, dim, dim> delta_sigma_mat;
    voigt_h::voigtToMatrix(delta_sigma, delta_sigma_mat);
    sigma_guess += delta_sigma_mat;

    projection_error = update_f(k_tmp, alpha_tmp);
    iterations++;
  }
}

/* -------------------------------------------------------------------------- */
/// Infinitesimal deformations
template <Int dim>
template <class Args>
inline void MaterialDruckerPrager<dim>::computeStressOnQuad(Args && args) {
  const auto & grad_u = args["grad_u"_n];
  const auto & previous_grad_u = args["previous_grad_u"_n];
  const auto & previous_sigma = args["previous_sigma"_n];
  const auto & sigma_th = args["sigma_th"_n];
  const auto & previous_sigma_th = args["previous_sigma_th"_n];

  Real delta_sigma_th = sigma_th - previous_sigma_th;

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;

  // Compute trial stress, sigma_tr
  Matrix<Real> sigma_tr(dim, dim);
  MaterialElastic<dim>::computeStressOnQuad(
      tuple::make_named_tuple("grad_u"_n = grad_delta_u, "sigma"_n = sigma_tr,
                              "sigma_th"_n = delta_sigma_th));
  sigma_tr += previous_sigma;

  bool initial_yielding = (this->computeYieldFunction(sigma_tr) > 0);

  // use closet point projection to compute the plastic multiplier and
  // gradient  and inealstic strain at the surface for the given trial stress
  // state
  Matrix<Real, dim, dim> delta_inelastic_strain;
  delta_inelastic_strain.zero();

  if (initial_yielding) {
    constexpr auto size = voigt_h::size;

    // plastic multiplier
    Real dp{0.};

    // gradient of yield surface in voigt notation
    Vector<Real, size> gradient_f;

    // inelastic strain in voigt notation
    Vector<Real, size> delta_inelastic_strain_voigt;

    // compute using closet-point projection
    this->computeGradientAndPlasticMultplier(sigma_tr, dp, gradient_f,
                                             delta_inelastic_strain_voigt);

    for (auto i : arange(dim, size)) {
      delta_inelastic_strain_voigt[i] /= 2.;
    }

    voigt_h::voigtToMatrix(delta_inelastic_strain_voigt,
                           delta_inelastic_strain);
  }

  // Compute the increment in inelastic strain
  MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(
      tuple::append(args, "delta_grad_u"_n = grad_delta_u,
                    "delta_inelastic_strain"_n = delta_inelastic_strain));
}

} // namespace akantu

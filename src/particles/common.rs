use super::Particle;

pub fn calculate_spin(particles: &mut [Particle]) {
    for (i, particle) in particles.iter_mut().enumerate() {
        assert_ne!(
            particle.moment_of_inertia, 0.,
            "Moment of inertia for particle {i} is zero!"
        );
        particle.spin.x = particle.angular_momentum.x / particle.moment_of_inertia;
        particle.spin.y = particle.angular_momentum.y / particle.moment_of_inertia;
        particle.spin.z = particle.angular_momentum.z / particle.moment_of_inertia;
        // norm needed for rotational flattening (torque and accelerations) and evolution
        particle.norm_spin_vector_2 =
            (particle.spin.x.powi(2)) + (particle.spin.y.powi(2)) + (particle.spin.z.powi(2));
    }
}

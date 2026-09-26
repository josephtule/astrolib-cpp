// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#include "core/body.hpp"

#include "core/ephemeris_provider.hpp"
#include "core/status.hpp"

#include <utility>

StatusCode set_celestial_ephemeris_providers(
    Celestial& celestial,
    BodyEphemerisProviders providers
) {
    StatusCode status;
    if (providers.translation != nullptr) {
        status = validate_ephemeris_provider(*providers.translation);
        if (status != StatusCode::ok) return status;
    }

    if (providers.orientation != nullptr) {
        status = validate_orientation_provider(*providers.orientation);
        if (status != StatusCode::ok) return status;
    }
    if (providers.orientation == nullptr
        && celestial.attitude_model == CelestialAttitudeModel::provider) {
        return StatusCode::invalid_input;
    }

    if (providers.translation != nullptr) celestial.propagate_tr = false;
    if (providers.orientation != nullptr) {
        celestial.attitude_model = CelestialAttitudeModel::provider;
        celestial.propagate_att = true;
    }
    celestial.ephemeris_providers = std::move(providers);

    return StatusCode::ok;
}

StatusCode instrument_suite_from_body(Body& body, InstrumentSuite*& out) {
    out = nullptr;
    switch (body.body_type) {
    case BodyType::station: {
        out = &static_cast<Station&>(body).instrument_suite;
    } break;
    case BodyType::satellite: {
        out = &static_cast<Satellite&>(body).instrument_suite;
    } break;
    case BodyType::celestial: [[fallthrough]];
    case BodyType::unknown: {
        return StatusCode::unsupported_type;
    } break;
    }

    return StatusCode::ok;
}

StatusCode instrument_suite_from_body(const Body& body, const InstrumentSuite*& out) {
    out = nullptr;
    switch (body.body_type) {
    case BodyType::station: {
        out = &static_cast<const Station&>(body).instrument_suite;
    } break;
    case BodyType::satellite: {
        out = &static_cast<const Satellite&>(body).instrument_suite;
    } break;
    case BodyType::celestial: [[fallthrough]];
    case BodyType::unknown: {
        return StatusCode::unsupported_type;
    } break;
    }

    return StatusCode::ok;
}

StatusCode mass_properties_from_body(Body& body, MassProperties*& out) {
    out = nullptr;
    switch (body.body_type) {
    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: return StatusCode::unsupported_type;
    case BodyType::satellite: {
        out = &static_cast<Satellite&>(body).mass_properties;
    } break;
    case BodyType::station: {
        out = &static_cast<Station&>(body).mass_properties;
    } break;
    }

    return StatusCode::ok;
}

StatusCode mass_properties_from_body(const Body& body, const MassProperties*& out) {
    out = nullptr;
    switch (body.body_type) {
    case BodyType::unknown: return StatusCode::unsupported_type;
    case BodyType::celestial: return StatusCode::unsupported_type;
    case BodyType::satellite: {
        out = &static_cast<const Satellite&>(body).mass_properties;
    } break;
    case BodyType::station: {
        out = &static_cast<const Station&>(body).mass_properties;
    } break;
    }

    return StatusCode::ok;
}

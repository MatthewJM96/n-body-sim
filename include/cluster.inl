#include <limits>
#include <random>

template <typename Precision, typename = typename std::enable_if<std::is_floating_point<Precision>::value>::type>
Precision cluster::member_distance_2(Member<Precision> lhs, Member<Precision> rhs) {
    return (lhs.x - rhs.x) * (lhs.x - rhs.x)
            + (lhs.y - rhs.y) * (lhs.y - rhs.y)
            + (lhs.z - rhs.z) * (lhs.z - rhs.z);
}

template <typename Precision>
void cluster::kpp(Member<Precision>* members, ui32 member_count, Member<Precision>* centroids, ui32 centroid_count) {
    /************
       Set up metadata for k++ algorithm.
                                ************/

    ui32 number_of_members_chosen = 1;
    bool* member_chosen_as_centroid = new bool[member_count];

    Precision* cumulative_distance_2s = new Precision[member_count];
    Precision  total_distance_2 = 0.0;

    std::default_random_engine generator;

    /************
       Make initial choice of a centroid.
                                ************/

    {
        std::uniform_int_distribution<ui32> distribution(0, member_count - 1);
        ui32 initial_choice = distribution(generator);

        centroids[0] = members[initial_choice];
        member_chosen_as_centroid[initial_choice] = true;
    }

    /************
       Perform k++ algorithm for each subsequent centroid.
                                                 ************/

    for (ui32 centroid_idx = 1; centroid_idx < centroid_count; ++centroid_idx) {
        //
        // For each member of the dataset not so far chosen as a centroid, determine
        // the minimum distance to a chosen centroid and select one of those members
        // to be the next centroid with probability proportional to distance^2 from
        // nearest centroid.
        //

        // Calculate distances for each member not so far chosen as a centroid.
        for (ui32 member_idx = 0; member_idx < member_count; ++member_idx) {
            // Skip members already chosen as a centroid.
            if (member_chosen_as_centroid[member_idx]) {
                cumulative_distance_2s[member_idx] = cumulative_distance_2s[member_idx - 1];
                continue;
            }

            // Calculate distance to nearest chosen centroid.
            Precision minimum_distance_2_to_chosen_centroids = std::numeric_limits<Precision>::max();
            for (ui32 chosen_centroid_idx = 0; chosen_centroid_idx < centroid_idx; ++chosen_centroid_idx) {
                Precision distance_2_to_chosen_centroid = member_distance_2(members[member_idx], centroids[chosen_centroid_idx]);

                if (distance_2_to_chosen_centroid < minimum_distance_2_to_chosen_centroids) {
                    minimum_distance_2_to_chosen_centroids = distance_2_to_chosen_centroid;
                }
            }

            // Place the calculated distance into metadata.
            if (member_idx == 0) {
                cumulative_distance_2s[member_idx] = minimum_distance_2_to_chosen_centroids;
            } else {
                cumulative_distance_2s[member_idx] = cumulative_distance_2s[member_idx - 1] + minimum_distance_2_to_chosen_centroids;
            }
            total_distance_2 += minimum_distance_2_to_chosen_centroids;
        }

        // Select a member to be next chosen centroid with probability proportional to
        // distance^2 to nearest existing centroid.
        std::uniform_real_distribution<Precision> distribution(0.0, total_distance_2);
        Precision choice = distribution(generator);

        for (ui32 member_idx = 0; member_idx < member_count; ++member_idx) {
            if (cumulative_distance_2s[member_idx] > choice) {
                centroids[centroid_idx] = members[member_idx];
                member_chosen_as_centroid[member_idx] = true;
            }
        }
    }
}

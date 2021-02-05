#include <algorithm>
#include <limits>
#include <random>

template <typename Precision>
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

template <typename Precision>
ui32 cluster::impl::nearest_centroid(const Member<Precision>& member, Member<Precision>* centroids, ui32 centroid_count) {
    ui32      nearest_centroid = 0;
    Precision nearest_centroid_distance_2 = std::numeric_limits<Precision>::max();

    for (ui32 centroid_idx = 0; centroid_idx < centroid_count; ++centroid_idx) {
        Precision centroid_distance_2 = member_distance_2(member, centroids[centroid_idx]);

        if (centroid_distance_2 < nearest_centroid_distance_2) {
            nearest_centroid = centroid_idx;
            nearest_centroid_distance_2 = centroid_distance_2;
        }
    }

    return nearest_centroid;
}

// TODO(Matthew): Eventually we want to be able to provide an initial cluster set and only do copies on entries that change from one cluster to another after so many iterations.
template <typename Precision>
void cluster::k_means(Member<Precision>* members, ui32 member_count, Member<Precision>* centroids, ui32 centroid_count, KMeansOptions options, OUT Cluster<Precision>*& clusters) {
    /************
       Set up metadata for k-means algorithm.
                                    ************/

    Member<Precision>* new_centroids = new Member<Precision>[centroid_count];
    ui32* member_centroid_map = new ui32[member_count];
    ui32* members_count_in_centroid = new ui32[member_count];

    /************
       Perform k-means algorithm.
                        ************/

    ui32 iterations = 0;
    ui32 changes_in_iteration = 0;
    do {
        // Complete if max iterations has been reached.
        if (++iterations > options.max_iterations) break;

        // Zero pertinent data.
        changes_in_iteration = 0;
        std::fill_n(new_centroids, centroid_count, Member<Precision>{0.0, 0.0, 0.0});
        std::fill_n(members_count_in_centroid, member_count, 0);

        // Iterate each member of the population on which the clusters are being built.
        // For each member, determine which centroid it is nearest to and add its position
        // to a new centroid which will then take its position as the average position of
        // the associated members.
        for (ui32 member_idx = 0; member_idx < member_count; ++member_idx) {
            ui32 nearest_centroid_idx = impl::nearest_centroid(members[member_idx], centroids, centroid_count);

            if (member_centroid_map[member_idx] != nearest_centroid_idx) {
                ++changes_in_iteration;
                member_centroid_map[member_idx] = nearest_centroid_idx;
            }

            new_centroids[nearest_centroid_idx].x += members[member_idx].x;
            new_centroids[nearest_centroid_idx].y += members[member_idx].y;
            new_centroids[nearest_centroid_idx].z += members[member_idx].z;

            ++members_count_in_centroid[nearest_centroid_idx];
        }

        // Using total members associated with each centroid, calculate the new
        // centroid for that group by taking the average of their positions.
        for (ui32 centroid_idx = 0; centroid_idx < centroid_count; ++centroid_idx) {
            new_centroids[centroid_idx].x /= members_count_in_centroid[centroid_idx];
            new_centroids[centroid_idx].y /= members_count_in_centroid[centroid_idx];
            new_centroids[centroid_idx].z /= members_count_in_centroid[centroid_idx];
        }
    } while (changes_in_iteration > options.acceptable_changes_per_iteration);
}

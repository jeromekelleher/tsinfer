#include "tsinfer.h"
#include "err.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

static void
ancestor_store_check_state(ancestor_store_t *self)
{
    int ret;
    site_id_t l, start, end, *focal_sites;
    size_t j, k;
    size_t total_segments = 0;
    size_t max_site_segments = 0;
    size_t num_older_ancestors, num_focal_sites, num_epoch_ancestors;
    ancestor_id_t *epoch_ancestors = malloc(self->num_ancestors * sizeof(ancestor_id_t));
    allele_t *a = malloc(self->num_sites * sizeof(allele_t));
    assert(a != NULL);
    assert(epoch_ancestors != NULL);

    for (l = 0; l < self->num_sites; l++) {
        total_segments += self->sites[l].num_segments;
        if (self->sites[l].num_segments > max_site_segments) {
            max_site_segments = self->sites[l].num_segments;
        }
    }
    assert(total_segments == self->total_segments);
    assert(max_site_segments == self->max_num_site_segments);
    for (j = 0; j < self->num_ancestors; j++) {
        ret = ancestor_store_get_ancestor(self, j, a, &start, &end,
                &num_older_ancestors, &num_focal_sites, &focal_sites);

        assert(ret == 0);
        /* if (j > 0) { */
        /*     assert(a[focal] == 1); */
        /*     /1* assert(start <= focal); *1/ */
        /*     assert(focal < end); */
        /* } */
        assert(end <= self->num_sites);
        assert(start < end);
        for (l = 0; l < self->num_sites; l++) {
            if (l < start || l >= end) {
                assert(a[l] == -1);
            } else {
                assert(a[l] != -1);
            }
        }
    }
    for (j = 1; j < self->num_epochs; j++) {
        ret = ancestor_store_get_epoch_ancestors(self, j,
                epoch_ancestors, &num_epoch_ancestors);
        assert(ret == 0);
        assert(num_epoch_ancestors > 0);
        for (k = 0; k < num_epoch_ancestors; k++) {
            assert(self->ancestors.age[epoch_ancestors[0]] ==
                    self->ancestors.age[epoch_ancestors[k]]);
        }
    }

    free(epoch_ancestors);
    free(a);
}

int
ancestor_store_print_state(ancestor_store_t *self, FILE *out)
{
    site_id_t l;
    site_state_t *site;
    size_t j, k;

    fprintf(out, "Ancestor store\n");
    fprintf(out, "num_sites = %d\n", (int) self->num_sites);
    fprintf(out, "num_ancestors = %d\n", (int) self->num_ancestors);
    fprintf(out, "total_segments  = %d\n", (int) self->total_segments);
    fprintf(out, "max_num_site_segments = %d\n", (int) self->max_num_site_segments);
    fprintf(out, "total_memory = %d\n", (int) self->total_memory);
    for (l = 0; l < self->num_sites; l++) {
        site = &self->sites[l];
        printf("%d\t%.3f\t[%d]:: ", (int) l, site->position, (int) site->num_segments);
        for (j = 0; j < site->num_segments; j++) {
            printf("(%d, %d)", site->start[j], site->end[j]);
        }
        printf("\n");
    }
    fprintf(out, "ancestors = \n");
    fprintf(out, "id\tage\tnum_older_ancestors\tnum_focal_sites\tfocal_sites\n");
    for (j = 0; j < self->num_ancestors; j++) {
        fprintf(out, "%d\t%d\t%d\t%d\t", (int) j,
                (int) self->ancestors.age[j],
                (int) self->ancestors.num_older_ancestors[j],
                (int) self->ancestors.num_focal_sites[j]);
        for (k = 0; k < self->ancestors.num_focal_sites[j]; k++) {
            fprintf(out, "%d", self->ancestors.focal_sites[j][k]);
            if (k < self->ancestors.num_focal_sites[j] - 1) {
                fprintf(out, ",");
            }
        }
        fprintf(out, "\n");
    }
    fprintf(out, "epochs = \n");
    fprintf(out, "id\tfirst_ancestor\tnum_ancestors\n");
    for (j = 0; j < self->num_epochs; j++) {
        fprintf(out, "%d\t%d\t%d\n", (int) j, (int) self->epochs.first_ancestor[j],
                (int) self->epochs.num_ancestors[j]);
    }
    ancestor_store_check_state(self);
    return 0;
}

int
ancestor_store_alloc(ancestor_store_t *self,
        size_t num_sites, double *position,
        size_t num_ancestors, uint32_t *ancestor_age,
        size_t num_focal_sites, ancestor_id_t *focal_site_ancestor, site_id_t *focal_site,
        size_t num_segments, site_id_t *site, ancestor_id_t *start, ancestor_id_t *end)
{
    int ret = 0;
    site_id_t j, l, site_start, site_end;
    size_t k, num_site_segments;
    uint32_t current_age, num_older_ancestors;
    size_t current_epoch;
    ancestor_id_t seg_num_ancestors, ancestor_id;

    /* TODO error checking */
    assert(num_sites > 0);
    assert(num_ancestors > 0);
    assert(num_focal_sites > 0);
    assert(num_focal_sites <= num_sites);
    assert(num_segments > 0);

    memset(self, 0, sizeof(ancestor_store_t));
    self->num_sites = num_sites;
    self->num_ancestors = num_ancestors;
    self->sites = calloc(num_sites, sizeof(site_state_t));
    self->ancestors.num_older_ancestors = malloc(num_ancestors * sizeof(uint32_t));
    self->ancestors.num_focal_sites = calloc(num_ancestors, sizeof(uint32_t));
    self->ancestors.focal_sites = malloc(num_ancestors * sizeof(site_id_t *));
    self->ancestors.age = malloc(num_ancestors * sizeof(uint32_t));
    self->ancestors.focal_sites_mem = malloc(num_focal_sites * sizeof(site_id_t));
    if (self->sites == NULL
            || self->ancestors.num_older_ancestors == NULL
            || self->ancestors.focal_sites == NULL
            || self->ancestors.num_focal_sites == NULL
            || self->ancestors.age == NULL
            || self->ancestors.focal_sites_mem == NULL) {
        ret = TSI_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(self->ancestors.age, ancestor_age, num_ancestors * sizeof(uint32_t));
    memcpy(self->ancestors.focal_sites_mem,
            focal_site, num_focal_sites * sizeof(site_id_t));
    self->ancestors.num_older_ancestors[0] = 0;
    self->ancestors.num_focal_sites[0] = 0;
    self->ancestors.focal_sites[0] = NULL;
    self->ancestors.age[0] = UINT32_MAX;
    ancestor_id = 0;
    for (k = 0; k < num_focal_sites; k++) {
        if (ancestor_id != focal_site_ancestor[k]) {
            ancestor_id++;
            /* TODO input checking here */
            assert(focal_site_ancestor[k] == ancestor_id);
            self->ancestors.focal_sites[ancestor_id] = self->ancestors.focal_sites_mem + k;
        }
        self->ancestors.num_focal_sites[ancestor_id]++;
    }

    site_start = 0;
    site_end = 0;
    self->max_num_site_segments = 0;
    seg_num_ancestors = 0;
    for (l = 0; l < self->num_sites; l++) {
        if (l > 0) {
            // TODO raise an error here.
            assert(position[l] > position[l - 1]);
        }
        self->sites[l].position = position[l];
        if (site_end < num_segments) {
            assert(site[site_start] >= l);
            assert(site[site_end] >= l);
            while (site_end < num_segments && site[site_end] == l) {
                site_end++;
            }
            num_site_segments = site_end - site_start;
            if (num_site_segments > self->max_num_site_segments) {
                self->max_num_site_segments = num_site_segments;
            }
            self->total_memory += num_site_segments * (2 * sizeof(ancestor_id_t) + sizeof(allele_t));
            self->sites[l].start = malloc(num_site_segments * sizeof(ancestor_id_t));
            self->sites[l].end = malloc(num_site_segments * sizeof(ancestor_id_t));
            if (self->sites[l].start == NULL || self->sites[l].end == NULL) {
                ret = TSI_ERR_NO_MEMORY;
                goto out;
            }
            k = 0;
            for (j = site_start; j < site_end; j++) {
                assert(site[j] == l);
                self->sites[l].start[k] = start[j];
                self->sites[l].end[k] = end[j];
                self->sites[l].num_segments++;
                self->total_segments++;
                if (end[j] > seg_num_ancestors) {
                    seg_num_ancestors = end[j];
                }
                k++;
            }
            site_start = site_end;
        }
    }

    /* Work out the number of epochs */
    self->num_epochs = 1;
    current_age = 0;
    for (j = 0; j < self->num_ancestors; j++) {
        if (self->ancestors.age[j] != current_age) {
            self->num_epochs++;
            current_age = self->ancestors.age[j];
        }
    }
    self->epochs.first_ancestor = calloc(self->num_epochs, sizeof(ancestor_id_t));
    self->epochs.num_ancestors = calloc(self->num_epochs, sizeof(size_t));
    if (self->epochs.first_ancestor == NULL || self->epochs.num_ancestors == NULL) {
        ret = TSI_ERR_NO_MEMORY;
        goto out;
    }
    /* Work out the number of older ancestors and assign to epochs */
    assert(self->num_ancestors > 1);
    current_epoch = self->num_epochs - 1;
    current_age = self->ancestors.age[1] + 1;
    self->epochs.first_ancestor[current_epoch] = 0;
    self->epochs.num_ancestors[current_epoch] = 1;
    num_older_ancestors = 0;
    for (j = 1; j < self->num_ancestors; j++) {
        if (self->ancestors.age[j] < current_age) {
            num_older_ancestors = j;
            current_age = self->ancestors.age[j];
            current_epoch--;
            self->epochs.first_ancestor[current_epoch] = j;
        }
        self->ancestors.num_older_ancestors[j] = num_older_ancestors;
        self->epochs.num_ancestors[current_epoch]++;
    }
    // TODO error checking.
    assert(self->total_segments == num_segments);
    assert(seg_num_ancestors == (ancestor_id_t) num_ancestors);
out:
    return ret;
}

int
ancestor_store_free(ancestor_store_t *self)
{
    site_id_t l;

    for (l = 0; l < self->num_sites; l++) {
        tsi_safe_free(self->sites[l].start);
        tsi_safe_free(self->sites[l].end);
    }
    tsi_safe_free(self->sites);
    tsi_safe_free(self->ancestors.num_older_ancestors);
    tsi_safe_free(self->ancestors.focal_sites);
    tsi_safe_free(self->ancestors.num_focal_sites);
    tsi_safe_free(self->ancestors.age);
    tsi_safe_free(self->ancestors.focal_sites_mem);
    tsi_safe_free(self->epochs.first_ancestor);
    tsi_safe_free(self->epochs.num_ancestors);
    return 0;
}

/*
 Returns the state of the specified ancestor at the specified site.
*/
int
ancestor_store_get_state(ancestor_store_t *self, site_id_t site_id,
        ancestor_id_t ancestor_id, allele_t *state)
{
    int ret = 0;
    site_state_t *site = &self->sites[site_id];
    size_t j = 0;

    j = 0;
    while (j < site->num_segments && site->end[j] <= ancestor_id) {
        j++;
    }
    *state = 0;
    if (j < site->num_segments &&
        site->start[j] <= ancestor_id && ancestor_id < site->end[j]) {
        *state = 1;
    }
    return ret;
}

int
ancestor_store_get_ancestor(ancestor_store_t *self, ancestor_id_t ancestor_id,
        allele_t *ancestor, site_id_t *start_site, site_id_t *end_site,
        size_t *num_older_ancestors,
        size_t *num_focal_sites, site_id_t **focal_sites)
{
    int ret = 0;
    site_id_t l, start;
    bool started = false;

    assert(ancestor_id < (ancestor_id_t) self->num_ancestors);
    memset(ancestor, 0xff, self->num_sites * sizeof(allele_t));
    start = 0;
    for (l = 0; l < self->num_sites; l++) {
        ret = ancestor_store_get_state(self, l, ancestor_id, ancestor + l);
        if (ret != 0) {
            goto out;
        }
        if (ancestor[l] != -1 && ! started) {
            start = l;
            started = true;
        }
        if (ancestor[l] == -1 && started) {
            break;
        }
    }
    *start_site = start;
    *focal_sites = self->ancestors.focal_sites[ancestor_id];
    *num_focal_sites = self->ancestors.num_focal_sites[ancestor_id];
    *num_older_ancestors = self->ancestors.num_older_ancestors[ancestor_id];
    *end_site = l;
out:
    return ret;
}


int
ancestor_store_get_epoch_ancestors(ancestor_store_t *self, int epoch,
        ancestor_id_t *epoch_ancestors, size_t *num_epoch_ancestors)
{
    int ret = 0;
    size_t j;

    assert(epoch > 0);
    assert(epoch < (int) self->num_epochs);
    for (j = 0; j < self->epochs.num_ancestors[epoch]; j++) {
        epoch_ancestors[j] = self->epochs.first_ancestor[epoch] + (ancestor_id_t) j;
    }
    *num_epoch_ancestors = self->epochs.num_ancestors[epoch];
    return ret;
}

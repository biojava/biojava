package org.biojava.nbio.structure.io.mmcif.model;

/**
 * Created by Matt on 10/31/2015.
 */
public class StructSiteGen extends AbstractBean {
    String id;
    String site_id;
    String auth_asym_id;
    String auth_atom_id;
    String auth_comp_id;
    String auth_seq_id;
    String label_alt_id;
    String label_asym_id;
    String label_atom_id;
    String label_comp_id;
    String label_seq_id;
    String details;
    String pdbx_auth_ins_code;
    String pdbx_num_res;
    String symmetry;

    public StructSiteGen() {
        super();
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getSite_id() {
        return site_id;
    }

    public void setSite_id(String site_id) {
        this.site_id = site_id;
    }

    public String getAuth_asym_id() {
        return auth_asym_id;
    }

    public void setAuth_asym_id(String auth_asym_id) {
        this.auth_asym_id = auth_asym_id;
    }

    public String getAuth_atom_id() {
        return auth_atom_id;
    }

    public void setAuth_atom_id(String auth_atom_id) {
        this.auth_atom_id = auth_atom_id;
    }

    public String getAuth_comp_id() {
        return auth_comp_id;
    }

    public void setAuth_comp_id(String auth_comp_id) {
        this.auth_comp_id = auth_comp_id;
    }

    public String getAuth_seq_id() {
        return auth_seq_id;
    }

    public void setAuth_seq_id(String auth_seq_id) {
        this.auth_seq_id = auth_seq_id;
    }

    public String getLabel_alt_id() {
        return label_alt_id;
    }

    public void setLabel_alt_id(String label_alt_id) {
        this.label_alt_id = label_alt_id;
    }

    public String getLabel_asym_id() {
        return label_asym_id;
    }

    public void setLabel_asym_id(String label_asym_id) {
        this.label_asym_id = label_asym_id;
    }

    public String getLabel_atom_id() {
        return label_atom_id;
    }

    public void setLabel_atom_id(String label_atom_id) {
        this.label_atom_id = label_atom_id;
    }

    public String getLabel_comp_id() {
        return label_comp_id;
    }

    public void setLabel_comp_id(String label_comp_id) {
        this.label_comp_id = label_comp_id;
    }

    public String getLabel_seq_id() {
        return label_seq_id;
    }

    public void setLabel_seq_id(String label_seq_id) {
        this.label_seq_id = label_seq_id;
    }

    public String getDetails() {
        return details;
    }

    public void setDetails(String details) {
        this.details = details;
    }

    public String getPdbx_auth_ins_code() {
        return pdbx_auth_ins_code;
    }

    public void setPdbx_auth_ins_code(String pdbx_auth_ins_code) {
        this.pdbx_auth_ins_code = pdbx_auth_ins_code;
    }

    public String getPdbx_num_res() {
        return pdbx_num_res;
    }

    public void setPdbx_num_res(String pdbx_num_res) {
        this.pdbx_num_res = pdbx_num_res;
    }

    public String getSymmetry() {
        return symmetry;
    }

    public void setSymmetry(String symmetry) {
        this.symmetry = symmetry;
    }
}

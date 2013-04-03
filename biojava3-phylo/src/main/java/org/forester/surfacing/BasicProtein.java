// $Id: BasicProtein.java,v 1.7 2009/10/26 23:29:40 cmzmasek Exp $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.surfacing;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

public class BasicProtein implements Protein {

    private final ProteinId    _id;
    private final Species      _species;
    private String             _name;
    private String             _desc;
    private String             _accession;
    private final List<Domain> _protein_domains;

    public BasicProtein( final String id_str, final String species_str ) {
        _id = new ProteinId( id_str );
        _species = new BasicSpecies( species_str );
        _protein_domains = new ArrayList<Domain>();
        init();
    }

    public void addProteinDomain( final Domain protein_domain ) {
        getProteinDomains().add( protein_domain );
    }

    @Override
    /**
     * If in_nc_order is set to true, this returns true only and only if
     * the order in List 'domains' and this protein (as determined by the start positions
     * of the domains of this proteins, _not_ by their index) are the same
     * (interspersing, 'other', domains in this are ignored). 
     * If in_nc_order is set to false, this returns true only and only if
     * this contains all domains listed in 'domains' (order and count do not matter).
     * 
     * @param domains a list of domain ids in a certain order.
     * @param in_nc_order to consider order
     * @return
     */
    public boolean contains( final List<DomainId> query_domain_ids, final boolean in_nc_order ) {
        if ( !in_nc_order ) {
            for( final DomainId query_domain_id : query_domain_ids ) {
                if ( !getProteinDomainIds().contains( query_domain_id ) ) {
                    return false;
                }
            }
            return true;
        }
        else {
            int current_start_position = -1;
            I: for( final DomainId query_domain_id : query_domain_ids ) {
                if ( getProteinDomainIds().contains( query_domain_id ) ) {
                    final List<Domain> found_domains = getProteinDomains( query_domain_id );
                    final SortedSet<Integer> ordered_start_positions = new TreeSet<Integer>();
                    for( final Domain found_domain : found_domains ) {
                        ordered_start_positions.add( found_domain.getFrom() );
                    }
                    for( final int start_position : ordered_start_positions ) {
                        if ( start_position > current_start_position ) {
                            current_start_position = start_position;
                            continue I;
                        }
                    }
                    return false;
                }
                else {
                    return false;
                }
            }
            return true;
        }
    }

    @Override
    public String getAccession() {
        return _accession;
    }

    @Override
    public String getDescription() {
        return _desc;
    }

    @Override
    public String getName() {
        return _name;
    }

    public int getNumberOfProteinDomains() {
        return getProteinDomains().size();
    }

    public Domain getProteinDomain( final int index ) {
        return _protein_domains.get( index );
    }

    public int getProteinDomainCount( final DomainId domain_id ) {
        return getProteinDomains( domain_id ).size();
    }

    private List<DomainId> getProteinDomainIds() {
        final List<DomainId> ids = new ArrayList<DomainId>( getProteinDomains().size() );
        for( final Domain domain : getProteinDomains() ) {
            ids.add( domain.getDomainId() );
        }
        return ids;
    }

    public List<Domain> getProteinDomains() {
        return _protein_domains;
    }

    public List<Domain> getProteinDomains( final DomainId domain_id ) {
        final List<Domain> domains = new ArrayList<Domain>();
        for( final Domain domain : getProteinDomains() ) {
            if ( domain.getDomainId().equals( domain_id ) ) {
                domains.add( domain );
            }
        }
        return domains;
    }

    public ProteinId getProteinId() {
        return _id;
    }

    public Species getSpecies() {
        return _species;
    }

    private void init() {
        _desc = "";
        _accession = "";
        _name = "";
    }

    public void setAccession( final String accession ) {
        _accession = accession;
    }

    public void setDescription( final String description ) {
        _desc = description;
    }

    public void setName( final String name ) {
        _name = name;
    }
}

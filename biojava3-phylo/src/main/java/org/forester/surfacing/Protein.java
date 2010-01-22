// $Id: Protein.java,v 1.6 2008/08/19 21:34:23 cmzmasek Exp $
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

import java.util.List;

public interface Protein {

    public void addProteinDomain( final Domain protein_domain );

    /**
     * If in_nc_order is set to true, this should return true only and only if
     * the order in List 'domains' and this protein (as determined by the start positions
     * of the domains of this proteins, _not_ by their index) are the same
     * (interspersing, 'other', domains in this are ignored). 
     * If in_nc_order is set to false, this should return true only and only if
     * this contains all domains listed in 'domains' (order and count do not matter).
     * 
     * @param domains a list of domain ids in a certain order.
     * @param in_nc_order to consider order
     * @return
     */
    public boolean contains( final List<DomainId> domains, final boolean in_nc_order );

    public String getAccession();

    public String getDescription();

    public String getName();

    public int getNumberOfProteinDomains();

    public Domain getProteinDomain( final int index );

    public int getProteinDomainCount( final DomainId domain_id );

    public List<Domain> getProteinDomains();

    public List<Domain> getProteinDomains( final DomainId domain_id );

    public ProteinId getProteinId();

    public Species getSpecies();
}
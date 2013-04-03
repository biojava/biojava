// $Id: PfamToGoMapping.java,v 1.6 2009/01/13 19:49:29 cmzmasek Exp $
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

package org.forester.go;

import org.forester.surfacing.DomainId;

public class PfamToGoMapping implements Mapping {

    private final DomainId _pfam_domain_id;
    private final GoId     _go_id;

    public PfamToGoMapping( final DomainId pfam_domain_id, final GoId go_id ) {
        _pfam_domain_id = pfam_domain_id;
        _go_id = go_id;
    }

    @Override
    public int compareTo( final Mapping m ) {
        if ( this == m ) {
            return 0;
        }
        return getKey().compareTo( ( DomainId ) m.getKey() );
    }

    /**
     * Based on key and value.
     * 
     * 
     */
    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check pfam to go mapping equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check dpfam to go mapping equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return getKey().equals( ( ( PfamToGoMapping ) o ).getKey() )
                    && getValue().equals( ( ( PfamToGoMapping ) o ).getValue() );
        }
    }

    @Override
    public DomainId getKey() {
        return _pfam_domain_id;
    }

    @Override
    public GoId getValue() {
        return _go_id;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getKey().toString() );
        sb.append( " > " );
        sb.append( getValue().toString() );
        return sb.toString();
    }
}

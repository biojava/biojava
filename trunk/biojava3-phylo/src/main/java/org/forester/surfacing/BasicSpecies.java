// $Id: BasicSpecies.java,v 1.5 2008/12/13 06:08:59 cmzmasek Exp $
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

import org.forester.util.ForesterUtil;

public class BasicSpecies implements Species {

    final private String _species_id;

    public BasicSpecies( final String species_id ) {
        if ( ForesterUtil.isEmpty( species_id ) ) {
            throw new IllegalArgumentException( "attempt to create new species from empty or null string" );
        }
        _species_id = species_id.trim();
    }

    @Override
    public int compareTo( final Species species ) {
        if ( this == species ) {
            return 0;
        }
        return getSpeciesId().toLowerCase().compareTo( species.getSpeciesId().toLowerCase() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return getSpeciesId().equals( ( ( Species ) o ).getSpeciesId() );
        }
    }

    /* (non-Javadoc)
     * @see org.forester.surfacing.Species#getSpeciesId()
     */
    public String getSpeciesId() {
        return _species_id;
    }

    @Override
    public int hashCode() {
        return getSpeciesId().hashCode();
    }

    @Override
    public String toString() {
        return getSpeciesId();
    }
}

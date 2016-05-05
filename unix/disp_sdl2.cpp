//******************************************************************************
///
/// @file unix/disp_sdl2.cpp
///
/// SDL2 (Simple direct media layer) based render display system.
///
/// @author Christoph Hormann <chris_hormann@gmx.de>
///
/// @copyright
/// @parblock
///
/// Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
/// Copyright 1991-2016 Persistence of Vision Raytracer Pty. Ltd.
///
/// POV-Ray is free software: you can redistribute it and/or modify
/// it under the terms of the GNU Affero General Public License as
/// published by the Free Software Foundation, either version 3 of the
/// License, or (at your option) any later version.
///
/// POV-Ray is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU Affero General Public License for more details.
///
/// You should have received a copy of the GNU Affero General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.
///
/// ----------------------------------------------------------------------------
///
/// POV-Ray is based on the popular DKB raytracer version 2.12.
/// DKBTrace was originally written by David K. Buck.
/// DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
///
/// @endparblock
///
//*******************************************************************************

#include "config.h"

#ifdef HAVE_LIBSDL2

#include "disp_sdl2.h"

#include <algorithm>

// this must be the last file included
#include "syspovdebug.h"


namespace pov_frontend
{
using namespace vfe;
using namespace vfePlatform;

extern shared_ptr<Display> gDisplay;

const UnixOptionsProcessor::Option_Info UnixSDL2Display::Options[] =
{
    // command line/povray.conf/environment options of this display mode can be added here
    // section name, option name, default, has_param, command line parameter, environment variable name, help text
    UnixOptionsProcessor::Option_Info("display", "scaled", "on", false, "", "POV_DISPLAY_SCALED", "scale render view to fit screen"),
    UnixOptionsProcessor::Option_Info("", "", "", false, "", "", "") // has to be last
};

bool UnixSDL2Display::Register(vfeUnixSession *session)
{
    session->GetUnixOptions()->Register(Options);
    // Initialize SDL
    if ( SDL_Init(SDL_INIT_VIDEO) != 0 )
    {
        fprintf(stderr, "Couldn't initialize SDL: %s.\n", SDL_GetError());
        return false;
    }

    atexit(SDL_Quit);
    // TODO: correct display detection
    return true;
}

UnixSDL2Display::UnixSDL2Display(unsigned int w, unsigned int h, GammaCurvePtr gamma, vfeSession *session, bool visible) :
    UnixDisplay(w, h, gamma, session, visible)
{
    m_window = NULL;
    m_renderer = NULL;
    m_surface = NULL;
    m_texture = NULL;
    m_thread = GetThreadId();
    m_state = F_SHOULD_DESTROY;
}

UnixSDL2Display::~UnixSDL2Display()
{
    Close();
}

void UnixSDL2Display::Initialise()
{
    if (m_VisibleOnCreation)
        Show();
}

void UnixSDL2Display::Hide()
{
}

bool UnixSDL2Display::TakeOver(UnixDisplay *display)
{
    UnixSDL2Display *p = dynamic_cast<UnixSDL2Display *>(display);
    if (p == NULL)
        return false;
    if ((GetWidth() != p->GetWidth()) || (GetHeight() != p->GetHeight()))
        return false;

    // Don't dead-lock.  (This really shouldn't be locked yet, anyways, since we -- presumably --
    // don't have anything yet!)
    boost::mutex::scoped_lock local_lock(control_mutex, boost::try_to_lock);
    if (!local_lock)
    {
        return false; // Ouch.  This really shouldn't happen, but...
    }
    boost::mutex::scoped_lock remote_lock(p->control_mutex);

    m_state = p->m_state;
    p->m_state &= ~F_SHOULD_DESTROY;

    m_thread = p->m_thread;
    m_texture = p->m_texture;
    m_renderer = p->m_renderer;
    m_window = p->m_window;
    m_surface = p->m_surface;
    m_update_rect = p->m_update_rect;
    m_PxCnt = p->m_PxCnt;

    return true;
}

void UnixSDL2Display::Close()
{
    boost::mutex::scoped_lock local_lock(control_mutex);

    if (m_state & F_SHOULD_DESTROY)
    {
        if (m_state & F_HAVE_SURFACE)
        {
            SDL_FreeSurface(m_surface);
        }
        if (m_state & F_HAVE_WINDOW)
        {
            SDL_DestroyRenderer(m_renderer);
            SDL_DestroyWindow(m_window);
        }
    }
    m_surface = NULL;
    m_renderer = NULL;
    m_texture = NULL;
    m_window = NULL;
    m_state = 0;

}

void UnixSDL2Display::Live()
{
    boost::mutex::scoped_lock local_lock(control_mutex);

    if ((m_state & F_HAVE_SURFACE) && !(m_state & F_HAVE_WINDOW))
    {
        m_thread = GetThreadId();

        Uint32 window_flags = SDL_WINDOW_ALLOW_HIGHDPI;
        int width = 0;
        int height = 0;

        vfeUnixSession *UxSession = dynamic_cast<vfeUnixSession *>(m_Session);

        m_display_scaled = UxSession->GetUnixOptions()->isOptionSet("display", "scaled");

        if (m_display_scaled)
            // determine maximum display area (wrong and ugly)
        {
            SDL_DisplayMode dm;
            if (SDL_GetDesktopDisplayMode(0, &dm) == 0)
            {
                {
                    width = min(dm.w - 10, (int)GetWidth());
                    height = min(dm.h - 80, (int)GetHeight());
                }

                // calculate display area
                float AspectRatio = float(width)/float(height);
                float AspectRatio_Full = float(GetWidth())/float(GetHeight());
                if (AspectRatio > AspectRatio_Full)
                    width = int(AspectRatio_Full*float(height));
                else if (AspectRatio != AspectRatio_Full)
                    height = int(float(width)/AspectRatio_Full);
            }
            else
            {
                fprintf(stderr, "Unable to determine display mode parameters: %s\n", SDL_GetError());
                return;
            }
        }
        else
        {
            window_flags |= SDL_WINDOW_FULLSCREEN_DESKTOP;
        }

        // Initialize the display
        SDL_CreateWindowAndRenderer(width, height, window_flags, &m_window, &m_renderer);

        if ( m_window == NULL )
        {
            if (!m_display_scaled)
            {
                fprintf(stderr, "Couldn't create fullscreen window: %s\n", SDL_GetError());
            }
            else
            {
                fprintf(stderr, "Couldn't create %dx%d window: %s\n", width, height, SDL_GetError());
            }
            return;
        }

        if ( m_renderer == NULL )
        {
            fprintf(stderr, "Couldn't create renderer: %s\n", SDL_GetError());
            return;
        }
        SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "best");

        if (SDL_RenderSetLogicalSize(m_renderer, GetWidth(), GetHeight()) != 0)
        {
            fprintf(stderr, "Couldn't set logical size of renderer to %dx%d: %s\n", GetWidth(), GetHeight(), SDL_GetError());
            return;
        }

        // Create texture
        m_texture = SDL_CreateTexture(m_renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, GetWidth(), GetHeight());

        if (m_texture == NULL)
        {
            fprintf(stderr, "Couldn't create %dx%d texture: %s\n", GetWidth(), GetHeight(), SDL_GetError());
            return;
        }

        m_state |= F_HAVE_WINDOW;
        SetCaption(false);
    }
    assert(m_thread == GetThreadId());
}



void UnixSDL2Display::SetCaption(bool paused)
{
    if (!(m_state & F_HAVE_WINDOW))
        return;

    boost::format f;
    if (m_display_scaled)
        f = boost::format(PACKAGE_NAME " " VERSION_BASE " SDL display (scaled)%s")
            % (paused ? " [paused]" : "");
    else
        f = boost::format(PACKAGE_NAME " " VERSION_BASE " SDL display%s")
            % (paused ? " [paused]" : "");
    // FIXME: SDL_WM_SetCaption() causes locks on some distros, see http://bugs.povray.org/23
    // FIXME: SDL_WM_SetCaption(f.str().c_str(), PACKAGE_NAME);
}

void UnixSDL2Display::Show()
{
    if (gDisplay.get() != this)
        gDisplay = m_Session->GetDisplay();

    boost::mutex::scoped_lock local_lock(control_mutex);

    if (!(m_state & F_HAVE_SURFACE))
    {

        // Create surface
        m_surface = SDL_CreateRGBSurface(0, GetWidth(), GetHeight(), 32, 0, 0, 0, 0);

        if (m_surface == NULL)
        {
            fprintf(stderr, "Couldn't create %dx%d surface: %s\n", GetWidth(), GetHeight(), SDL_GetError());
            return;
        }

        m_update_rect.x = 0;
        m_update_rect.y = 0;
        m_update_rect.w = GetWidth();
        m_update_rect.h = GetHeight();

        m_state |= F_HAVE_SURFACE;
        m_PxCnt = UpdateInterval;

    }
}

inline void UnixSDL2Display::SetPixel(unsigned int x, unsigned int y, const RGBA8& colour)
{
    *(Uint32 *)((Uint8 *)m_surface->pixels + y * m_surface->pitch + x * 4) = SDL_MapRGBA(m_surface->format, colour.red, colour.green, colour.blue, colour.alpha);
}

void UnixSDL2Display::UpdateCoord(unsigned int x, unsigned int y)
{
    unsigned int rx2 = m_update_rect.x + m_update_rect.w;
    unsigned int ry2 = m_update_rect.y + m_update_rect.h;
    m_update_rect.x = min((unsigned int)m_update_rect.x, x);
    m_update_rect.y = min((unsigned int)m_update_rect.y, y);
    rx2 = max(rx2, x);
    ry2 = max(ry2, y);
    m_update_rect.w = rx2 - m_update_rect.x;
    m_update_rect.h = ry2 - m_update_rect.y;
}

void UnixSDL2Display::UpdateCoord(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
    unsigned int rx2 = m_update_rect.x + m_update_rect.w;
    unsigned int ry2 = m_update_rect.y + m_update_rect.h;
    m_update_rect.x = min((unsigned int)m_update_rect.x, x1);
    m_update_rect.y = min((unsigned int)m_update_rect.y, y1);
    rx2 = max(rx2, x2);
    ry2 = max(ry2, y2);
    m_update_rect.w = rx2 - m_update_rect.x;
    m_update_rect.h = ry2 - m_update_rect.y;
}

void UnixSDL2Display::DrawPixel(unsigned int x, unsigned int y, const RGBA8& colour)
{
    if (!(m_state & F_HAVE_SURFACE) || x >= GetWidth() || y >= GetHeight())
        return;

    SetPixel(x, y, colour);
    UpdateCoord(x, y);

    m_PxCnt++;

}

void UnixSDL2Display::DrawRectangleFrame(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8& colour)
{
    if (!(m_state & F_HAVE_SURFACE))
        return;

    int ix1 = min(x1, GetWidth()-1);
    int ix2 = min(x2, GetWidth()-1);
    int iy1 = min(y1, GetHeight()-1);
    int iy2 = min(y2, GetHeight()-1);

    for(unsigned int x = ix1; x <= ix2; x++)
    {
        SetPixel(x, iy1, colour);
        SetPixel(x, iy2, colour);
    }

    for(unsigned int y = iy1; y <= iy2; y++)
    {
        SetPixel(ix1, y, colour);
        SetPixel(ix2, y, colour);
    }
    UpdateCoord(ix1, iy1, ix2, iy2);

    m_PxCnt = UpdateInterval;
}

void UnixSDL2Display::DrawFilledRectangle(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8& colour)
{
    if (!(m_state & F_HAVE_SURFACE))
        return;

    unsigned int ix1 = min(x1, GetWidth()-1);
    unsigned int ix2 = min(x2, GetWidth()-1);
    unsigned int iy1 = min(y1, GetHeight()-1);
    unsigned int iy2 = min(y2, GetHeight()-1);

    UpdateCoord(ix1, iy1, ix2, iy2);

    SDL_Rect tempRect;
    tempRect.x = ix1;
    tempRect.y = iy1;
    tempRect.w = ix2 - ix1 + 1;
    tempRect.h = iy2 - iy1 + 1;
    SDL_FillRect(m_surface, &tempRect, SDL_MapRGBA(m_surface->format, colour.red, colour.green, colour.blue, colour.alpha));

    m_PxCnt = UpdateInterval;
}

void UnixSDL2Display::DrawPixelBlock(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8 *colour)
{
    if (!(m_state & F_HAVE_SURFACE))
        return;

    unsigned int ix1 = min(x1, GetWidth()-1);
    unsigned int ix2 = min(x2, GetWidth()-1);
    unsigned int iy1 = min(y1, GetHeight()-1);
    unsigned int iy2 = min(y2, GetHeight()-1);

    for(unsigned int y = y1, i = 0; y <= iy2; y++)
        for(unsigned int x = ix1; x <= ix2; x++, i++)
            SetPixel(x, y, colour[i]);
    UpdateCoord(ix1, iy1, ix2, iy2);

    m_PxCnt = UpdateInterval;
}

void UnixSDL2Display::Clear()
{
    m_update_rect.x = 0;
    m_update_rect.y = 0;
    m_update_rect.w = GetWidth();
    m_update_rect.h = GetHeight();

    SDL_FillRect(m_surface, NULL, (Uint32)0);

    m_PxCnt = UpdateInterval;
}

void UnixSDL2Display::UpdateScreen(bool Force = false)
{
    if (!(m_state & F_HAVE_WINDOW))
        return;

    assert(m_thread == GetThreadId());
    if (Force || m_PxCnt >= UpdateInterval)
    {
        SDL_UpdateTexture(m_texture, &m_update_rect, m_surface->pixels, m_surface->pitch);
        SDL_RenderClear(m_renderer);
        SDL_RenderCopy(m_renderer, m_texture, NULL, NULL);
        SDL_RenderPresent(m_renderer);
        m_PxCnt = 0;
    }
}

void UnixSDL2Display::PauseWhenDoneNotifyStart()
{
    fprintf(stderr, "Press a key or click the display to continue...");
    SetCaption(true);
}

void UnixSDL2Display::PauseWhenDoneNotifyEnd()
{
    SetCaption(false);
    fprintf(stderr, "\n\n");
}

bool UnixSDL2Display::PauseWhenDoneResumeIsRequested()
{
    if (SDL_WasInit(SDL_INIT_EVENTS) == 0)
        return true;

    SDL_Event event;
    bool do_quit = false;

    if (SDL_PollEvent(&event))
    {
        switch (event.type)
        {
        case SDL_KEYDOWN:
            if ( event.key.keysym.sym == SDLK_q || event.key.keysym.sym == SDLK_RETURN || event.key.keysym.sym == SDLK_KP_ENTER )
                do_quit = true;
            break;
        case SDL_MOUSEBUTTONDOWN:
            do_quit = true;
            break;
        }
    }

    return do_quit;
}

bool UnixSDL2Display::HandleEvents()
{
    if (SDL_WasInit(SDL_INIT_EVENTS) == 0)
        return false;

    SDL_Event event;
    bool do_quit = false;

    while (SDL_PollEvent(&event))
    {
        switch (event.type)
        {
        case SDL_KEYDOWN:
            if ( event.key.keysym.sym == SDLK_q )
                do_quit = true;
            else if ( event.key.keysym.sym == SDLK_p )
            {
                if (!m_Session->IsPausable())
                    break;

                if (m_Session->Paused())
                {
                    if (m_Session->Resume())
                        SetCaption(false);
                }
                else
                {
                    if (m_Session->Pause())
                        SetCaption(true);
                }
            }
            break;
        case SDL_QUIT:
            do_quit = true;
            break;
        }
        if (do_quit)
            break;
    }

    return do_quit;
}

}

#endif /* HAVE_LIBSDL2 */

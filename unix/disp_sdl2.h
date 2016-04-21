/*******************************************************************************
 * disp_sdl2.h
 *
 * Written by Christoph Hormann <chris_hormann@gmx.de>
 *
 * SDL2 (Simple direct media layer) based render display system
 *
 * ---------------------------------------------------------------------------
 * Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
 * Copyright 1991-2013 Persistence of Vision Raytracer Pty. Ltd.
 *
 * POV-Ray is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * POV-Ray is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------------
 * POV-Ray is based on the popular DKB raytracer version 2.12.
 * DKBTrace was originally written by David K. Buck.
 * DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
 * ---------------------------------------------------------------------------
 * $File: //depot/public/povray/3.x/unix/disp_sdl.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifdef HAVE_LIBSDL2

#ifndef _DISP_SDL2_H
#define _DISP_SDL2_H

#include "vfe.h"
#include "unixoptions.h"
#include "disp.h"

#include <boost/thread.hpp>
#include <SDL2/SDL.h>

namespace pov_frontend
{
	using namespace vfe;
	using namespace vfePlatform;

	class UnixSDL2Display : public UnixDisplay
	{
		public:
			static const UnixOptionsProcessor::Option_Info Options[];
			static bool Register(vfeUnixSession *session);

			UnixSDL2Display(unsigned int w, unsigned int h, GammaCurvePtr gamma, vfeSession *session, bool visible);
			virtual ~UnixSDL2Display();
			void Initialise();
			void Close();
			void Show();
			void Hide();
			bool TakeOver(UnixDisplay *display);
			void DrawPixel(unsigned int x, unsigned int y, const RGBA8& colour);
			void DrawRectangleFrame(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8& colour);
			void DrawFilledRectangle(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8& colour);
			void DrawPixelBlock(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8 *colour);
			void Clear();
			bool HandleEvents();
			void Live();
			void UpdateScreen(bool Force);
			void PauseWhenDoneNotifyStart();
			bool PauseWhenDoneResumeIsRequested();
			void PauseWhenDoneNotifyEnd();

		protected:
			/// Number of Pixels before the display is updated
			static const unsigned int UpdateInterval = 100;

			/// State flags
			static const unsigned int F_HAVE_SURFACE = 0x01;
			static const unsigned int F_HAVE_WINDOW = 0x02;
			static const unsigned int F_SHOULD_DESTROY = 0x04;

			void SetCaption(bool paused);

			/// Sets the color of a pixel in a non-scaled image.
			inline void SetPixel(unsigned int x, unsigned int y, const RGBA8& colour);
			/// Makes a pixel coordinate part of the update rectangle
			void UpdateCoord(unsigned int x, unsigned int y);
			/// Makes a rectangle part of the update rectangle
			void UpdateCoord(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);

			unsigned int m_state;
			bool m_display_scaled;
			/// for update interval
			unsigned int m_PxCnt;
			SDL_Window *m_window;
			SDL_Renderer *m_renderer;
			SDL_Rect m_update_rect;
			SDL_Surface *m_surface;
			SDL_Texture *m_texture;
			POVMS_Sys_Thread_Type m_thread;
			boost::mutex control_mutex;
	};
}

#endif /* _DISP_SDL2_H */

#endif /* HAVE_LIBSDL2 */

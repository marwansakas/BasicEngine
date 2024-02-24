#pragma once
#ifndef TEXTURE_H
#define TEXTURE_H

#include <string>

class Texture
{
public:
	Texture(const std::string& fileName);
	Texture(int width, int height,unsigned char *data);
	void Bind(int slot);
	inline int GetSlot(){return m_texture;}
	 ~Texture();
protected:
private:
	Texture(const Texture& texture) {}
	void operator=(const Texture& texture) {}
	unsigned int m_texture;
};

#endif

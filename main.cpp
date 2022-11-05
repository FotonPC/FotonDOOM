//  FotonPC's school projects
//  october 2022
//
//





#include <SFML\\Graphics.hpp>
#include <SFML\\Audio.hpp>
#include "engine.h"
#include <vector>
#include <windows.h>
#include <algorithm>
#include <cctype>

sf::Uint8 convert_to_uint8(int x) {
    return sf::Uint8(min(255, max(0, x)));
};

void load_engine(Engine3D &engine, int ry, double ray_step_koef, int rx, int n_r_sp, int n_f_sp, int n_clls) {

    engine.init(2, 2, 90, "resources\\images\\bg_menu" + std::to_string(ry) + ".png", ray_step_koef, rx, ry, n_r_sp, n_f_sp, n_clls);
};


int main()

{   
    bool is_recommended_use_shaders = True;
    // Get GPU ID 
    DISPLAY_DEVICEA dd;
    dd.cb = sizeof(DISPLAY_DEVICEA);
    EnumDisplayDevicesA(NULL, 0, &dd, EDD_GET_DEVICE_INTERFACE_NAME);
    string gpu_name = string(dd.DeviceString);
    
    std::transform(gpu_name.begin(), gpu_name.end(), gpu_name.begin(),
        [](unsigned char c) { return std::tolower(c); });
    cout << gpu_name << endl;
    if (gpu_name.find(string("hd")) != std::string::npos) {
        is_recommended_use_shaders = False;
        cout << "Using shaders not recomended! Bad GPU" << endl;
    }
    int RESOLUTION_X = 512;
    int RESOLUTION_Y = 256;
    int rx = RESOLUTION_X;
    int ry = RESOLUTION_Y;
    bool is_use_shaders = True;
    
    double ray_step_koef;
    cout << "begin init" << endl;
    std::string line;
    std::ifstream in("settings.txt");
    int str_len, k2;
    float u_blur = 10;
    float u_viewnetka = 0.6;
    double light_force = 30;
    int screen_resolution_d_unstable;
    k2 = 0;
    if (in.is_open())
    {
        cout << "File with MAP is loaded" << endl;
        while (getline(in, line))
        {
            if (k2 == 0) {
                if (line == "low") {
                    ray_step_koef = 100;
                }
                if (line == "medium") {
                    ray_step_koef = 20;
                }
                if (line == "high") {
                    ray_step_koef = 10;
                }
                if (line == "ultra") {
                    ray_step_koef = 5;
                }
                if (line == "epic") {
                    ray_step_koef = 2;
                }
                if (line == "ultra epic") {
                    ray_step_koef = 1;
                }
            }
            if (k2 == 1) {
                if (line == "low") {
                    ry = 32;
                }
                if (line == "medium") {
                    ry = 64;
                }
                if (line == "high") {
                    ry = 128;
                }
                if (line == "ultra") {
                    ry = 256;
                }
                if (line == "epic") {
                    ry = 512;
                }
                if (line == "ultra epic") {
                    ry = 1024;
                }
            }
            if (k2 == 2) {
                screen_resolution_d_unstable = atoi(line.c_str());
            }
            if (k2 == 3) {
                if (line == "used") {
                    is_use_shaders = True;
                }
                if (line == "notused") {
                    is_use_shaders = False;
            }
            }
            if (k2 == 4) {
                u_blur = atof(line.c_str());
            }
            if (k2 == 5) {
                u_viewnetka = atof(line.c_str());
            }
            if (k2 == 6) {
                light_force = atof(line.c_str());
            }
            k2++;
        }
    }
    in.close();
    rx = ry * 2;
    RESOLUTION_X = rx;
    RESOLUTION_Y = ry;
    double real_width, real_height;
    real_width = rx;
    real_height = ry;
    
    sf::RenderWindow window(sf::VideoMode(RESOLUTION_X, RESOLUTION_Y), "AgentF");

    int n_r_sp = 2;
    int n_f_sp = 10;
    int n_clls = 1;

    sf::Font font;
    font.loadFromFile("resources\\sansation.ttf");
    sf::Text text("", font, 20);
    text.setStyle(sf::Text::Bold | sf::Text::Underlined);
    text.setString("Loading engine");
    text.setPosition(0, 0);
    window.draw(text);
    window.display();
    Engine3D engine;
    /*thread Thread_(&load_engine, &engine, ry, ray_step_koef, rx, n_r_sp, n_f_sp, n_clls);
    while (engine.not_load) {
        window.draw(text);
        window.display();
    }
    Thread_.join();*/
    engine.init(2, 2, 90, "resources\\images\\bg_menu" + std::to_string(ry) + ".png", ray_step_koef, rx, ry, n_r_sp, n_f_sp, n_clls);
    engine.light_force = light_force;
    engine.is_use_shaders = is_use_shaders;
    sf::Texture texture;
    sf::Texture h_map_texture;
    sf::RenderTexture texture2, h_map_render;
    texture2.create(RESOLUTION_X, RESOLUTION_Y);
    h_map_render.create(rx, ry);
    h_map_texture = h_map_render.getTexture();
    texture = texture2.getTexture();
    sf::Uint8* pixels = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];

    sf::Sprite sprite(texture);
    engine.page = 0;
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Engine loaded. Loading level");
    window.draw(text);
    window.display();
    engine.load_level("flat1.level");
    //engine.set_map("resources\\maps\\flat1.map");
    engine.load_texture(0, "resources\\images\\textures\\parquet"+ std::to_string(ry)+".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\parquet" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(1, "resources\\images\\textures\\brick" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\brick" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(11, "resources\\images\\textures\\up" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\up" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(2, "resources\\images\\textures\\wood" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\wood" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(3, "resources\\images\\textures\\marble" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\marble" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(15, "resources\\images\\textures\\plika_" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\plika_" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(12, "resources\\images\\sprites\\door_1\\1_" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\sprites\\door_1\\1_" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(13, "resources\\images\\sprites\\door_1\\2_" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\sprites\\door_1\\2_" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(14, "resources\\images\\sprites\\door_1\\1_2_" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\sprites\\door_1\\1_2_" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(4, "resources\\images\\textures\\shukat_" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\shukat_" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(9, "resources\\images\\textures\\lift_" + std::to_string(ry) + ".png");
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Loaded resources\\images\\textures\\lift_" + std::to_string(ry) + ".png");
    window.draw(text);
    window.display();
    engine.load_texture(16, "resources\\images\\textures\\shukat_lamp_" + std::to_string(ry) + ".png");
    //engine.load_texture(4, "resources\\images\\textures\\stone_old" + std::to_string(ry) + ".png");
    //engine.load_texture(5, "resources\\images\\textures\\door" + std::to_string(ry) + ".png");
    //engine.load_texture(10, "resources\\images\\textures\\fog" + std::to_string(ry) + ".png");
    sf::Image icon;
    if (!icon.loadFromFile("resources\\logo.png"))
    {
        return 1;
    }
    window.setIcon(512, 512, icon.getPixelsPtr());
    window.clear();
    pixels = engine.alternative_stage(texture, sprite, 255);
    texture.update(pixels);
    window.draw(sprite);
    text.setString("Level & Textures loaded. Loading music and finish initizialing");
    window.draw(text);
    window.display();

    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    sf::Vertex point;
    
    window.setFramerateLimit(120);
    engine.page = 0;
    sf::Mouse mouse;
    mouse = sf::Mouse();
    sf::SoundBuffer buffer1;
    buffer1.loadFromFile("resources\\audio\\parquet.flac");
    sf::Sound sound1;
    sound1.setBuffer(buffer1);
    sound1.setVolume(5);
    sf::SoundBuffer buffer2;
    buffer2.loadFromFile("resources\\audio\\door_open.flac");
    sf::Sound sound2;
    sound2.setBuffer(buffer2);
    sound2.setVolume(50);
    sf::SoundBuffer buffer3;
    buffer3.loadFromFile("resources\\audio\\door_close.ogg");
    sf::Sound sound3;
    sound3.setBuffer(buffer3);
    sound3.setVolume(50);
    bool is_step;
    int step_counter = 0;
    int door_close_counter = 0;
    int door_open_counter = 0;

    sf::Music music_butt1;
    music_butt1.openFromFile("resources\\audio\\butt1.ogg");

    sf::Music music_menu;
    music_menu.openFromFile("resources\\audio\\menu_v4.ogg");
    music_menu.play();
    music_menu.setVolume(50);
    bool is_butt1_play;
    is_butt1_play = False;
    int alpha_after_menu = 0;
    


    sf::Clock clock;
    int i, j, k;
    int iteration_counter;
    iteration_counter = 0;
    double fps;
    double a;
        a = engine.sprites_objects[0].x_size;
        engine.sprites_objects[0].x_size = engine.sprites_objects[0].y_size;
        engine.sprites_objects[0].y_size = a;
        engine.sprites_objects[0].orient_x = 1 - engine.sprites_objects[0].orient_x;

        if (engine.sprites_objects[0].y == 9.05) {
            engine.sprites_objects[0].y = 9.25;
            engine.sprites_objects[0].x = 4.95;
            door_open_counter = 1000;
        }
        else {

            engine.sprites_objects[0].y = 9.05;
            engine.sprites_objects[0].x = 4.75;
            door_close_counter = 1000;

        }

        engine.flat_sprites[0].load_texture(1, "resources\\images\\hands2_" + std::to_string(ry) + ".png");

        sf::RenderTexture outputTexture;
        sf::Shader shader;
        sf::Sprite outputTextureSpriteFlipped;
        sf::Sprite outputTextureSprite;

        if (is_use_shaders) {

            
            outputTexture.create(rx, ry);
            outputTextureSprite = sf::Sprite(outputTexture.getTexture());
            outputTextureSpriteFlipped = sf::Sprite(texture);
            outputTextureSpriteFlipped.setScale(1, -1);
            outputTextureSpriteFlipped.setPosition(0, ry);
            shader.loadFromFile("shader.frag", sf::Shader::Fragment);
            shader.setUniform("rx", rx);
            shader.setUniform("ry", ry);
            shader.setUniform("u_blur", u_blur);
            shader.setUniform("u_viewnetka", u_viewnetka);
        }
        
       
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            
            if (event.type == sf::Event::MouseMoved)
            {
                if (engine.page == 1) {
                    try {
                        engine.angle -= sf::Mouse::getPosition(window).x - (int)window.getSize().x / 2 + screen_resolution_d_unstable;
                        while (engine.angle >= 360) {
                            engine.angle -= 360.0;
                        }
                        //engine.angle = (long long)engine.angle % 360;
                        if (engine.angle < 0) {
                            engine.angle = to_double(360) + engine.angle;
                        }
                        //cout << sf::Mouse::getPosition(window).x - (int)window.getSize().x / 2 +  8 << endl;

                        mouse.setPosition(sf::Vector2i(window.getPosition().x + real_width / 2, window.getPosition().y + real_height / 2));
                    }
                    catch (...) {
                        cout << "Error in mouse move" << endl;
                    }
                    engine.flat_sprites[0].set_pos(engine.x + cos(engine.angle / 180 * PI) / 2, engine.y + sin(engine.angle / 180 * PI) / 2);
                    engine.flat_sprites[0].set_angle(engine.angle / 180 * PI + PI / 2);
                }
            }
            if (event.type == sf::Event::Resized) {
                real_width = event.size.width;
                real_height = event.size.height;
            }
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Escape) {
                    if (engine.page == 0) {
                        window.close();
                    }
                }
                if (event.key.code == sf::Keyboard::Space) {
                    cout << engine.x << " " << engine.y << " " << engine.angle << endl;
                    double a;
                    if (door_open_counter > 1.5 * fps && door_close_counter > fps / 2 && engine.x < 5 && engine.x > 4 && engine.y > 8.3 && engine.y < 10) {
                        a = engine.sprites_objects[0].x_size;
                        double x_size_past = engine.sprites_objects[0].x_size;
                        double y_size_past = engine.sprites_objects[0].y_size;
                        bool orient_x_past = engine.sprites_objects[0].orient_x;
                        double x_past = engine.sprites_objects[0].x;
                        double y_past = engine.sprites_objects[0].y;
                        engine.sprites_objects[0].x_size = engine.sprites_objects[0].y_size;
                        engine.sprites_objects[0].y_size = a;
                        engine.sprites_objects[0].orient_x = 1 - engine.sprites_objects[0].orient_x;

                        if (engine.sprites_objects[0].y == 9.05) {
                            engine.sprites_objects[0].y = 9.25;
                            engine.sprites_objects[0].x = 4.95;
                            sound2.play();
                            door_open_counter = 0;
                        }
                        else {

                            engine.sprites_objects[0].y = 9.05;
                            engine.sprites_objects[0].x = 4.75;
                            sound3.play();
                            door_close_counter = 0;
                        }
                        if (engine.sprites_objects[0].is_collision(engine.x, engine.y)) {
                            engine.sprites_objects[0].x_size = x_size_past;
                            engine.sprites_objects[0].y_size = y_size_past;
                            engine.sprites_objects[0].orient_x = orient_x_past;
                            engine.sprites_objects[0].y = y_past;
                            engine.sprites_objects[0].x = x_past;
                            sound2.stop();
                            sound3.stop();
                        }
                    }
                    
                    
                   
                }
                if (event.key.code == sf::Keyboard::F) {
                    if (engine.is_butt_1) {
                        if (not is_butt1_play) {
                            music_butt1.play();
                            is_butt1_play = true;
                        }
                        else {
                            music_butt1.pause();
                            is_butt1_play = False;
                        }
                    }
                }
                if (event.key.code == sf::Keyboard::Enter) {
                    if (engine.page == 0) {
                        engine.is_stoped_heights = false;
                        engine.page = 1;
                        alpha_after_menu = 1;
                        window.setMouseCursorVisible(false);
                        mouse.setPosition(sf::Vector2i(window.getPosition().x + real_width / 2, window.getPosition().y + real_height / 2));
                        window.setMouseCursorGrabbed(true);
                        music_menu.stop();
                        //engine.stage_heights();


                    }
                    else if (engine.page == 1) {
                        engine.is_stoped_heights = true;
                        engine.page = 0;
                        alpha_after_menu = 0;
                        window.setMouseCursorVisible(true);
                        window.setMouseCursorGrabbed(false);
                        music_butt1.stop();
                        music_menu.play();
                        engine.is_stoped_heights = true;
                    }
                }
            }
        }
        sf::Time time = clock.getElapsedTime();
        fps = 1.0f / time.asSeconds();
        try {
            if (fps > 1 && iteration_counter % to_int(fps) < fps / 5.0) {
                window.setTitle("AgentF: FPS: "+to_string(to_int(fps)));
            }
        }
        catch (...) {

        }
        
        clock.restart().asSeconds();

        
        if (engine.page == 1) {
            is_step = false;
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.angle += 1;
                //cout << 120 / fps << endl;
                //cout << engine.angle << endl;
                while (engine.angle >= 360) {
                    engine.angle -= 360.0;
                }
                //engine.angle = (long long)engine.angle % 360;
                if (engine.angle < 0) {
                    engine.angle = to_double(360) + engine.angle;
                }
                //cout << sf::Mouse::getPosition(window).x - (int)window.getSize().x / 2 +  8 << endl;
                

            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.angle -=1;
                while (engine.angle >= 360) {
                    engine.angle -= 360.0;
                }
                //engine.angle = (long long)engine.angle % 360;
                if (engine.angle < 0) {
                    engine.angle = to_double(360) + engine.angle;
                }
                //cout << sf::Mouse::getPosition(window).x - (int)window.getSize().x / 2 +  8 << endl;

                cout << engine.angle << endl;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.step_left(4/fps);

                if (step_counter > fps / 3.5) {
                    sound1.stop();
                    step_counter = 0;
                    sound1.play();
                    engine.flat_sprites[0].ctc = 1 - engine.flat_sprites[0].ctc;
                }
                is_step = true;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.step_right(4 / fps);
                if (step_counter > fps / 3.5) {
                    sound1.stop();
                    step_counter = 0;
                    sound1.play();
                    engine.flat_sprites[0].ctc = 1 - engine.flat_sprites[0].ctc;
                }
                is_step = true;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.step_forward(4 / fps);
                if (step_counter > fps / 3.5) {
                    sound1.stop();
                    step_counter = 0;
                    sound1.play();
                    engine.flat_sprites[0].ctc = 1 - engine.flat_sprites[0].ctc;
                }
                is_step = true;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.step_backward(4 / fps);
                if (step_counter > fps / 3.5) {
                    sound1.stop();
                    step_counter = 0;
                    sound1.play();
                    engine.flat_sprites[0].ctc = 1 - engine.flat_sprites[0].ctc;
                }
                is_step = true;
            }
            if (not is_step) {
                step_counter = fps;
            }
            while (engine.angle >= 360) {
                engine.angle -= 360.0;
            }
            if (engine.angle < 0) {
                engine.angle += 360.0;
            }
            begin = std::chrono::steady_clock::now();
            //engine.thread_stage_heights(0, rx, 1);
            //engine.stage_heights();
            end = std::chrono::steady_clock::now();
            elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
            //if (iteration_counter % 30 == 0) { std::cout << "Heights time: " << elapsed_ms.count() << " mcs\n"; }
            music_butt1.setVolume(1/(sqrt((engine.x-1.6) * (engine.x - 1.6) + (engine.y - 2.9) * (engine.y - 2.9))) * 20);
            step_counter++;
            door_open_counter++;
            door_close_counter++;
            
        }
        window.clear();
        begin = std::chrono::steady_clock::now();
        //engine.thread_stage_render(0, rx, 1, 255);
        //pixels = engine.pixels;
        
       //pixels = engine.stage_render(texture, sprite, alpha_after_menu);
        engine.flat_sprites[0].set_pos(engine.x + cos(engine.angle / 180 * PI) / 2, engine.y + sin(engine.angle / 180 * PI) / 2);
        engine.flat_sprites[0].set_angle(engine.angle / 180 * PI + PI / 2);
        engine.flat_sprites[9].set_angle(engine.angle / 180 * PI + PI / 2);
        engine.flat_sprites[9].set_diff((engine.x-engine.flat_sprites[9].x-0.1)/fps, (engine.y-engine.flat_sprites[9].y-0.1)/fps);
        //pixels = engine.alternative_stage(texture, sprite, alpha_after_menu);
        pixels = engine.alternative_stage(texture, sprite, alpha_after_menu);
        texture.update(pixels);
        if (is_use_shaders && engine.page == 1) {
            pixels = engine.return_h_map();
            h_map_texture.update(pixels);
        }
        end = std::chrono::steady_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
        if (iteration_counter % 30 == 0) { std::cout << "Render stage time: " << elapsed_ms.count() << " mcs\n"; }

        begin = std::chrono::steady_clock::now();
        if (is_use_shaders && engine.page==1) {
            
            shader.setUniform("u_resolution", sf::Glsl::Vec2{ window.getSize() });
            shader.setUniform("u_time", clock.getElapsedTime().asSeconds());
            shader.setUniform("u_sample", texture);
            shader.setUniform("u_h_map", h_map_texture);
            window.draw(sprite, &shader);
        }
        else {
            window.draw(sprite);
        }
        end = std::chrono::steady_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
        if (iteration_counter % 30 == 0) { std::cout << "Shader time: " << elapsed_ms.count() << " mcs\n"; }
        
        /*
#pragma omp parallel for

        for (i = 0; i < rx; i++) {
#pragma omp parallel for
            for (j = 0; j < ry; j++) {
                //cout << engine.render_matrix[i][j][0] << endl;
                pixels[j * rx * 4 + i * 4] = engine.render_matrix[i][j][0];
                pixels[j * rx * 4 + i * 4 + 1] = engine.render_matrix[i][j][1];
                pixels[j * rx * 4 + i * 4 + 2] = engine.render_matrix[i][j][2];
                pixels[j * rx * 4 + i * 4 + 3] = 255;
            }
        }
        texture.update(pixels);
        */
        

        /*
    #pragma omp parallel for 
        for (i = 0; i < rx; i++) {
        #pragma omp parallel for 
            for (j = 0; j < ry; j++) {
                //cout << engine.render_matrix[i][j][0] << endl;
                point.position = sf::Vector2f(i, j);
                point.color = sf::Color{ convert_to_uint8(engine.render_matrix[i][j][0]),convert_to_uint8(engine.render_matrix[i][j][1]),convert_to_uint8(engine.render_matrix[i][j][2]) };
                window.draw(&point, 1, sf::Points);
            }
        }*/
        
        window.display();
        iteration_counter++;
        engine.is_butt_1 = (1.5 < engine.x && 2 > engine.x) && (2.5 < engine.y && 3 > engine.y);
        alpha_after_menu = min(255, alpha_after_menu + 1);
    }

    return 0;
};

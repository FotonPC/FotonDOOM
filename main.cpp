#include <SFML\\Graphics.hpp>
#include "engine.h"
#include <vector>

sf::Uint8 convert_to_uint8(int x) {
    return sf::Uint8(min(255, max(0, x)));
};




int main()

{
    int RESOLUTION_X = 512;
    int RESOLUTION_Y = 256;
    int rx = RESOLUTION_X;
    int ry = RESOLUTION_Y;
    
    double ray_step_koef;
    cout << "begin init" << endl;
    std::string line;
    std::ifstream in("settings.txt");
    int str_len, k2;
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
    
    sf::RenderWindow window(sf::VideoMode(RESOLUTION_X, RESOLUTION_Y), "FotonDOOM on C++");
    
    Engine3D engine;
    engine.init(1.5, 1.5, 90, "resources\\images\\bg" + std::to_string(ry) + ".png", ray_step_koef, rx, ry);
    engine.set_map("level.map");
    engine.load_texture(0, "resources\\images\\textures\\parquet"+ std::to_string(ry)+".png");
    engine.load_texture(1, "resources\\images\\textures\\brick" + std::to_string(ry) + ".png");
    engine.load_texture(11, "resources\\images\\textures\\up" + std::to_string(ry) + ".png");
    engine.load_texture(2, "resources\\images\\textures\\wood" + std::to_string(ry) + ".png");
    engine.load_texture(3, "resources\\images\\textures\\marble" + std::to_string(ry) + ".png");
    //engine.load_texture(4, "resources\\images\\textures\\stone_old" + std::to_string(ry) + ".png");
    //engine.load_texture(5, "resources\\images\\textures\\door" + std::to_string(ry) + ".png");
    //engine.load_texture(10, "resources\\images\\textures\\fog" + std::to_string(ry) + ".png");

    

    cout << "File with MAP is loaded 3" << endl;
    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    sf::Vertex point;
    sf::Texture texture;
    texture.create(RESOLUTION_X, RESOLUTION_Y);
    sf::Uint8* pixels = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];

    sf::Sprite sprite(texture);
    window.setFramerateLimit(3000);
    engine.page = 0;
    sf::Mouse mouse;
    mouse = sf::Mouse();


    sf::Clock clock;
    int i, j, k;
    int iteration_counter;
    iteration_counter = 0;
    double fps;
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
                        engine.angle -= sf::Mouse::getPosition(window).x - (int)window.getSize().x / 2 + 8;
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
                }
            }
            if (event.type == sf::Event::Resized) {
                real_width = event.size.width;
                real_height = event.size.height;
            }
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Enter) {
                    if (engine.page == 0) {
                        engine.page = 1;
                        window.setMouseCursorVisible(false);
                        mouse.setPosition(sf::Vector2i(window.getPosition().x + real_width / 2, window.getPosition().y + real_height / 2));
                        window.setMouseCursorGrabbed(true);

                    }
                    else if (engine.page == 1) {
                        engine.page = 0;
                        window.setMouseCursorVisible(true);
                        window.setMouseCursorGrabbed(false);
                    }
                }
            }
        }
        sf::Time time = clock.getElapsedTime();
        if (iteration_counter % 60 == 0) {
            cout << 1.0f / time.asSeconds() << endl;
        }
        fps = 1.0f / time.asSeconds();
        clock.restart().asSeconds();

        window.clear();
        if (engine.page == 1) {
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
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.step_right(4 / fps);
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.step_forward(4 / fps);
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
            {
                // левая стрелка нажата - двигаем нашего персонажа
                engine.step_backward(4 / fps);
            }
            while (engine.angle >= 360) {
                engine.angle -= 360.0;
            }
            if (engine.angle < 0) {
                engine.angle += 360.0;
            }
            begin = std::chrono::steady_clock::now();
            engine.stage_heights();
            end = std::chrono::steady_clock::now();
            elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
            //if (iteration_counter % 30 == 0) { std::cout << "Heights time: " << elapsed_ms.count() << " mcs\n"; }
        }
        begin = std::chrono::steady_clock::now();
        pixels = engine.stage_render(texture, sprite);
        texture.update(pixels);
        window.draw(sprite);
        end = std::chrono::steady_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
        //if (iteration_counter % 30 == 0) { std::cout << "Render stage time: " << elapsed_ms.count() << " mcs\n"; }
        begin = std::chrono::steady_clock::now();
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
        end = std::chrono::steady_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
        //if (iteration_counter % 30 == 0) { std::cout << "Render time: " << elapsed_ms.count() << " mcs\n"; }
        window.display();
        iteration_counter++;
    }

    return 0;
};
